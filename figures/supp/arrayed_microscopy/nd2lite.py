"""A small read-only ND2 reader, stdlib + numpy only.

Why this exists: no ND2 reader is installed anywhere on this Mac (`nd2`,
`nd2reader`, `pims`, `aicsimageio`, Fiji, bioformats - none), and installing one
was not available when this figure was built. The LV2.34 plates are plain
uncompressed ND2, so the container is cheap to walk directly.

Format, as verified against the six LV2.34 files:

* Every chunk is `uint32 magic=0x0ABECEDA | uint32 name_len | uint64 data_len |
  name | data`.
* The last 40 bytes carry `ND2 CHUNK MAP SIGNATURE 0000001!` followed by a
  uint64 pointing at the chunk-map chunk.
* Inside the chunk map each entry is `name | uint64 offset | uint64 length`.
  **The name carries no terminator.** Offsets are 4096-aligned, so their low
  byte is always 0x00 and reads exactly like a trailing NUL - scan to that byte
  and then start the offset *at* it, not after it. Getting this wrong shifts
  every name by one character and every offset into nonsense.
* An `ImageDataSeq|N!` chunk is `uint64 lead | pixels`, with the components
  interleaved per pixel and a row stride of `uiWidthBytes`. `_selfcheck` asserts
  the length arithmetic, which is what catches a wrong lead or a compressed file.

Metadata is Nikon's CLx variant blob. Rather than implement that, the two
helpers below pull values out positionally: a scalar or a string stored
immediately after its UTF-16LE key name. That is enough for the width, height,
component count, pixel calibration and stage positions, and it fails loudly
(empty list) rather than quietly when a key is absent.
"""
import os
import struct

import numpy as np

MAGIC = 0x0ABECEDA
SIG = b"ND2 CHUNK MAP SIGNATURE 0000001!"
LEAD = 8          # uint64 ahead of the pixels in every ImageDataSeq chunk


class ND2:
    def __init__(self, path):
        self.path = str(path)
        self.fh = open(self.path, "rb")
        self.map = self._read_map()
        attrs = self.raw("ImageAttributesLV!")
        self.width = scalar_after(attrs, "uiWidth")[0]
        self.height = scalar_after(attrs, "uiHeight")[0]
        self.comp = scalar_after(attrs, "uiComp")[0]
        self.row_bytes = scalar_after(attrs, "uiWidthBytes")[0]
        self.bpc = scalar_after(attrs, "uiBpcInMemory")[0]
        self.um_per_px = scalar_after(
            self.raw("ImageCalibrationLV|0!"), "dCalibration", "<d")[0]

    # -- container ---------------------------------------------------------
    def _chunk(self, offset):
        self.fh.seek(offset)
        magic, nl, dl = struct.unpack("<IIQ", self.fh.read(16))
        if magic != MAGIC:
            raise ValueError(f"{self.path}: bad chunk magic {magic:#x} at {offset}")
        self.fh.read(nl)
        return offset + 16 + nl, dl

    def _read_map(self):
        self.fh.seek(-40, os.SEEK_END)
        tail = self.fh.read(40)
        i = tail.find(SIG)
        if i < 0:
            raise ValueError(f"{self.path}: no chunk-map signature")
        (moff,) = struct.unpack("<Q", tail[i + len(SIG):i + len(SIG) + 8])
        dstart, dlen = self._chunk(moff)
        self.fh.seek(dstart)
        data = self.fh.read(dlen)
        out, p = {}, 0
        while True:
            e = data.find(b"\x00", p)
            if e < 0 or e + 16 > len(data):
                break
            name = data[p:e].decode("ascii", "replace")
            if not name or name.startswith("ND2 CHUNK MAP SIGNATURE"):
                break
            out[name] = struct.unpack("<QQ", data[e:e + 16])
            p = e + 16
        return out

    def raw(self, key):
        off, _ = self.map[key]
        dstart, dlen = self._chunk(off)
        self.fh.seek(dstart)
        return self.fh.read(dlen)

    # -- pixels ------------------------------------------------------------
    def n_frames(self):
        return sum(1 for k in self.map if k.startswith("ImageDataSeq|"))

    def _pixels(self, frame):
        dstart, dlen = self._chunk(self.map[f"ImageDataSeq|{frame}!"][0])
        expect = self.height * self.row_bytes + LEAD
        if dlen != expect:
            raise ValueError(
                f"{self.path} frame {frame}: {dlen} bytes, expected {expect}; "
                "the file is not plain uncompressed ND2")
        return dstart + LEAD

    def crop(self, frame, size, y0=None, x0=None):
        """A `size`x`size`x`comp` uint16 window, centred unless told otherwise.

        Read as one contiguous band of full rows and sliced in memory: over SMB
        that is ~30x faster than seeking per row (11 MB/s against 0.35).
        """
        base = self._pixels(frame)
        y0 = (self.height - size) // 2 if y0 is None else y0
        x0 = (self.width - size) // 2 if x0 is None else x0
        y0 = max(0, min(y0, self.height - size))
        x0 = max(0, min(x0, self.width - size))
        self.fh.seek(base + y0 * self.row_bytes)
        band = self.fh.read(size * self.row_bytes)
        a = np.frombuffer(band, "<u2").reshape(size, self.width, self.comp)
        return a[:, x0:x0 + size, :].copy()

    def coarse(self, frame, step):
        """The whole frame subsampled every `step` px, as uint16.

        One seek per sampled row. Slower per byte than `crop` but it returns the
        entire well, which is the only way to see that a well has a cell-free
        centre rather than no cells.
        """
        base = self._pixels(frame)
        rows = []
        for y in range(0, self.height, step):
            self.fh.seek(base + y * self.row_bytes)
            row = np.frombuffer(self.fh.read(self.row_bytes), "<u2")
            rows.append(row.reshape(self.width, self.comp)[::step, :])
        return np.array(rows)

    # -- acquisition layout ------------------------------------------------
    def frame_wells(self):
        """Frame index -> well name, matched on recorded stage coordinates.

        Not by position in the point list: plate 6 stops after 20 of the 24
        listed points, so index order would silently mislabel 16 wells.
        """
        meta = self.raw("ImageMetadataLV!")
        names = string_after(meta, "dPosName")
        px = scalar_after(meta, "dPosX", "<d")
        py = scalar_after(meta, "dPosY", "<d")
        xs = struct.unpack(f"<{len(self.raw('CustomData|X!')) // 8}d",
                           self.raw("CustomData|X!"))
        ys = struct.unpack(f"<{len(self.raw('CustomData|Y!')) // 8}d",
                           self.raw("CustomData|Y!"))
        out = []
        for x, y in zip(xs, ys):
            d = [((a - x) ** 2 + (b - y) ** 2) ** 0.5 for a, b in zip(px, py)]
            j = int(np.argmin(d))
            if d[j] > 100.0:          # wells sit 19.6 mm apart; 100 um is slack
                raise ValueError(f"{self.path}: frame at {x},{y} matches no point")
            out.append(names[j])
        return out

    def channels(self):
        """Component order, from the per-channel exposure chunk names."""
        pre = "CustomData|Camera_ExposureTime1_"
        hits = []
        for k in self.map:
            if k.startswith(pre):
                rest = k[len(pre):].rstrip("!")
                idx, _, name = rest.partition("_")
                if idx.isdigit():
                    hits.append((int(idx), name))
        return [n for _, n in sorted(hits)]


# -- positional metadata helpers ------------------------------------------
def scalar_after(b, key, fmt="<I"):
    """Every scalar stored straight after a UTF-16LE `key` and its NUL."""
    pat = key.encode("utf-16-le") + b"\x00\x00"
    size = struct.calcsize(fmt)
    out, start = [], 0
    while True:
        i = b.find(pat, start)
        if i < 0:
            return out
        j = i + len(pat)
        if j + size <= len(b):
            out.append(struct.unpack(fmt, b[j:j + size])[0])
        start = i + 1


def string_after(b, key):
    """Every UTF-16LE string stored straight after a UTF-16LE `key` and its NUL."""
    pat = key.encode("utf-16-le") + b"\x00\x00"
    out, start = [], 0
    while True:
        i = b.find(pat, start)
        if i < 0:
            return out
        j = i + len(pat)
        e = b.find(b"\x00\x00", j)
        if e < 0:
            return out
        if (e - j) % 2:
            e += 1
        out.append(b[j:e].decode("utf-16-le", "replace"))
        start = i + 1


def _selfcheck(path):
    """Fails if the container walk, the frame arithmetic or the well match break."""
    n = ND2(path)
    assert n.row_bytes == n.width * n.comp * n.bpc // 8, "row stride disagrees"
    assert n.n_frames() == scalar_after(
        n.raw("ImageAttributesLV!"), "uiSequenceCount")[0], "frame count disagrees"
    for f in range(n.n_frames()):
        _, dlen = n._chunk(n.map[f"ImageDataSeq|{f}!"][0])
        assert dlen == n.height * n.row_bytes + LEAD, f"frame {f} length"
    wells = n.frame_wells()
    assert len(wells) == n.n_frames() == len(set(wells)), "well match is not 1:1"
    assert 0.1 < n.um_per_px < 10, "implausible calibration"
    print(f"ok  {os.path.basename(path)}  {n.width}x{n.height}x{n.comp} "
          f"{n.bpc}-bit  {n.n_frames()} frames  {n.um_per_px:.4f} um/px  "
          f"{n.channels()}  {wells[0]}..{wells[-1]}")


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        _selfcheck(p)
