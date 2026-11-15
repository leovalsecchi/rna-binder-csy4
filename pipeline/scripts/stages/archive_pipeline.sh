#!/usr/bin/env bash
# archive_pipeline.sh — pipeline-safe wrapper around the lab's tar.zst utils.
#
# Adds the three things the raw tarz/untarz lack for use inside the rna_binder
# pipeline (concurrent jobs, shared dirs, crash-safety):
#   1. flock   — pack/unpack a unit under an exclusive per-unit lock, so
#                --num-runs>1 jobs sharing one designs/ and parallel refolds
#                never race on the same unit.
#   2. atomic  — tarz writes to <unit>.tar.zst.tmp, then mv into place; the
#                source dir is deleted ONLY after the verified archive exists.
#   3. sidecar — a plain <unit>.tar.zst.list (file list) written at pack time,
#                so listing/counting never has to stream-decompress the archive
#                (tar.zst has no central directory; this is what zip gave free).
#
# A "unit" is a stage-output dir: e.g. <run>/designs, <run>/refold_af3.
# At rest a unit is the single file <unit>.tar.zst (+ .list, + .lock).
#
# Verbs:
#   pack_unit   <dir> [tarz opts]     dir -> dir.tar.zst (+ .list); deletes dir
#   unpack_unit <archive> [dest]      archive -> dest/<base>/ (default: in place)
#   ensure_unit <dir>                 make <dir>/ exist (unpack in place if needed)
#   extract_member <archive> <glob> [dest]   pull matching members (ChimeraX open)
#   count_unit  <dir|archive>         #files, from the .list sidecar (no decompress)
#   purge       <path>...             rm -rf disposables, guarded
#
# Requires: bash, flock, zstd, tar, realpath, mktemp, find  (all on Euler).
# Source the lab utils (tarz/untarz/lsz/xtarz) — overridable via $ARCHIVE_SH.

: "${ARCHIVE_SH:=/cluster/work/wicky/utils/files/archive.sh}"
if ! command -v tarz >/dev/null 2>&1; then
    if [[ -f "$ARCHIVE_SH" ]]; then
        # shellcheck source=/dev/null
        source "$ARCHIVE_SH"
    else
        echo "archive_pipeline: tarz not found and ARCHIVE_SH ($ARCHIVE_SH) missing" >&2
        return 1 2>/dev/null || exit 1
    fi
fi

# pack_unit <dir> [extra tarz opts, e.g. -l 19 -T 8]
# dir -> dir.tar.zst + dir.tar.zst.list, then rm -rf dir. flock'd, atomic,
# verified-before-delete (tarz runs zstd -t). No-op if already packed.
pack_unit() {
    local dir="${1%/}"; shift || true
    local archive="${dir}.tar.zst" tmp="${dir}.tar.zst.tmp"
    local sidecar="${dir}.tar.zst.list" lock="${dir}.lock"
    (
        flock -x 200 || { echo "pack_unit: could not lock $lock" >&2; exit 1; }
        if [[ ! -d "$dir" && -f "$archive" ]]; then exit 0; fi          # already packed
        if [[ ! -d "$dir" ]]; then echo "pack_unit: $dir absent, no archive" >&2; exit 1; fi
        # Member list from the live dir (free — no decompress) for the sidecar.
        ( cd "$(dirname "$dir")" && find "$(basename "$dir")" ! -type d | sort ) > "${sidecar}.tmp" \
            || { echo "pack_unit: listing failed for $dir" >&2; rm -f "${sidecar}.tmp"; exit 1; }
        rm -f "$tmp"
        # tarz parses leading opts, then archive ($tmp); -C makes member paths relative.
        tarz "$@" "$tmp" -C "$(dirname "$dir")" "$(basename "$dir")" \
            || { echo "pack_unit: tarz failed for $dir (kept originals)" >&2; rm -f "$tmp" "${sidecar}.tmp"; exit 1; }
        mv -f "$tmp" "$archive" && mv -f "${sidecar}.tmp" "$sidecar" \
            || { echo "pack_unit: publish failed for $archive" >&2; exit 1; }
        rm -rf "$dir"
        echo "pack_unit: $dir -> $archive ($(wc -l < "$sidecar") files)"
    ) 200>"$lock"
}

# unpack_unit <archive.tar.zst> [dest_dir]
# Extract into dest_dir (default: the archive's own dir) as <dest>/<base>/.
# flock'd, atomic (extract to tmp then mv), idempotent (no-op if already there).
unpack_unit() {
    local archive; archive="$(realpath -m "$1")"
    local dest="${2:-$(dirname "$archive")}"
    local base; base="$(basename "${archive%.tar.zst}")"
    local lock="${archive%.tar.zst}.lock"
    (
        flock -x 200 || { echo "unpack_unit: could not lock $lock" >&2; exit 1; }
        [[ -d "$dest/$base" ]] && exit 0                               # already expanded
        [[ -f "$archive" ]] || { echo "unpack_unit: no archive $archive" >&2; exit 1; }
        mkdir -p "$dest"
        local tmp; tmp="$(mktemp -d "$dest/.unpack.XXXXXX")" || exit 1
        ( cd "$tmp" && untarz "$archive" ) \
            || { echo "unpack_unit: extract failed for $archive" >&2; rm -rf "$tmp"; exit 1; }
        mv "$tmp/$base" "$dest/$base" && rm -rf "$tmp" \
            || { echo "unpack_unit: publish failed for $dest/$base" >&2; rm -rf "$tmp"; exit 1; }
        echo "unpack_unit: $archive -> $dest/$base"
    ) 200>"$lock"
}

# ensure_unit <dir> : guarantee <dir>/ exists, unpacking <dir>.tar.zst in place.
ensure_unit() {
    local dir="${1%/}"
    [[ -d "$dir" ]] && return 0
    unpack_unit "${dir}.tar.zst" "$(dirname "$dir")"
}

# extract_member <archive.tar.zst> <member-glob> [dest_dir]
# Pull specific members (e.g. one design for ChimeraX). NOTE: tar.zst is
# sequential, so this still stream-decompresses the whole archive — fine for a
# per-run designs archive given zstd's speed; cache to /tmp for repeat opens.
extract_member() {
    local archive; archive="$(realpath -m "$1")"
    local glob="$2" dest="${3:-.}"
    mkdir -p "$dest"
    ( cd "$dest" && xtarz "$archive" "$glob" )
}

# count_unit <dir|archive> : number of files, from the .list sidecar only.
count_unit() {
    local base="${1%/}"; base="${base%.tar.zst}"
    local sidecar="${base}.tar.zst.list"
    if [[ -f "$sidecar" ]]; then wc -l < "$sidecar"; else echo 0; fi
}

# purge <path>... : rm -rf inode-heavy disposables, guarded.
# Caller expands globs (e.g. `purge "$RUN"/mpnn_native_design*`); each arg must
# live under a campaigns/ tree and not be empty/root/$HOME.
purge() {
    local p
    for p in "$@"; do
        [[ -e "$p" ]] || continue
        case "$p" in ""|"/"|"$HOME"|"$HOME/") echo "purge: refusing '$p'" >&2; continue ;; esac
        [[ "$p" == *"/campaigns/"* ]] || { echo "purge: refusing non-campaign '$p'" >&2; continue; }
        rm -rf -- "$p" && echo "purge: removed $p"
    done
}

# ── per-uuid slice verbs (added 2026-06-09 for --num-runs>1 design packing) ──
# A "slice" is one run's share of a SHARED designs/ dir, identified by its uuid8
# token. Each of N concurrent runs packs only its own *_<token>_* files into
# <designs_dir>/<token>.tar.zst, so the shared summary dir ends up holding N
# slice archives (zero loose files) and each refold unpacks just its own slice.
# Race-safe: a slice op only ever touches *_<token>_* members + its own
# <token>.tar.zst* sidecars, so sibling runs never collide. tar|zstd (not tarz)
# to match the readers (archive_read.py / aggregate.py both `zstd -d | tar -x`).

# pack_slice <designs_dir> <token>
# Pack <designs_dir>/*_<token>_* (flat members) -> <designs_dir>/<token>.tar.zst
# (+ .list), then rm the now-archived loose files. flock'd, atomic, verified.
# No-op if already packed or no matching members.
pack_slice() {
    local dir="${1%/}" tok="$2"
    local arch="$dir/$tok.tar.zst" tmp="$dir/$tok.tar.zst.tmp"
    local sidecar="$dir/$tok.tar.zst.list" lock="$dir/$tok.tar.zst.lock"
    (
        flock -x 200 || { echo "pack_slice: could not lock $lock" >&2; exit 1; }
        if [[ -f "$arch" ]]; then
            # [Fix 3a] already published; finish any loose-member cleanup a crash
            # between the mv and the rm below may have left half-done (idempotent;
            # the sidecar lists exactly this slice's archived members).
            [[ -f "$sidecar" ]] && ( cd "$dir" && tr '\n' '\0' < "$sidecar" | xargs -0 rm -f )
            exit 0
        fi
        ( cd "$dir" && find . -maxdepth 1 -type f -name "*_${tok}_*" \
              ! -name '*.tar.zst' ! -name '*.tar.zst.list' -printf '%P\n' | sort ) > "$sidecar.tmp" \
            || { echo "pack_slice: listing failed in $dir" >&2; rm -f "$sidecar.tmp"; exit 1; }
        if [[ ! -s "$sidecar.tmp" ]]; then
            rm -f "$sidecar.tmp"; echo "pack_slice: no *_${tok}_* members in $dir"; exit 0
        fi
        rm -f "$tmp"
        ( cd "$dir" && tar -c --files-from "$sidecar.tmp" | zstd -q -19 -T0 -o "$tmp" ) \
            || { echo "pack_slice: tar|zstd failed for $tok (kept originals)" >&2; rm -f "$tmp" "$sidecar.tmp"; exit 1; }
        zstd -tq "$tmp" \
            || { echo "pack_slice: archive verify failed for $tmp" >&2; rm -f "$tmp" "$sidecar.tmp"; exit 1; }
        mv -f "$tmp" "$arch" && mv -f "$sidecar.tmp" "$sidecar" \
            || { echo "pack_slice: publish failed for $arch" >&2; exit 1; }
        ( cd "$dir" && tr '\n' '\0' < "$sidecar" | xargs -0 rm -f )  # drop the now-archived loose files
        echo "pack_slice: $dir [*_${tok}_*] -> $arch ($(wc -l < "$sidecar") files)"
    ) 200>"$lock"
}

# unpack_slice <slice.tar.zst> <dest_dir>
# Expand a per-uuid slice archive (flat members) into dest_dir. flock'd.
unpack_slice() {
    local arch; arch="$(realpath -m "$1")"
    local dest="${2:?unpack_slice needs a dest dir}"
    local lock="${arch}.lock"
    (
        flock -x 200 || { echo "unpack_slice: could not lock $lock" >&2; exit 1; }
        [[ -f "$arch" ]] || { echo "unpack_slice: no archive $arch" >&2; exit 1; }
        mkdir -p "$dest"
        zstd -dq -c "$arch" | ( cd "$dest" && tar -x ) \
            || { echo "unpack_slice: extract failed for $arch" >&2; exit 1; }
        echo "unpack_slice: $arch -> $dest"
    ) 200>"$lock"
}
