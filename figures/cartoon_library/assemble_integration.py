"""Lay every re-themed piece out as one figure.

Takes the files `retheme_integration.py` writes and places them: the three rows,
a cross on each att site, and the +Bxb1 arrow between the landing pad and the
integrated locus. Every piece in the set is used, so nothing can be left behind.

Each file keeps its own .cls- names, so they are namespaced per piece on the way
in; otherwise the last style block would win and recolour everything.

    python code/cartoon_library/assemble_integration.py <dir-with-the-pieces>
"""
import re, sys
from pathlib import Path
S = Path(sys.argv[1]); d = S if (S / "Promoterless-donor-plasmid.svg").exists() else S / "zip_out"

def part(name):
    t = (d / name).read_text()
    vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', t).group(1).split()]
    style = re.search(r'<style>(.*?)</style>', t, re.S).group(1)
    body = re.sub(r'^.*?</defs>', '', t, flags=re.S)
    return vb, style, body[:body.rfind("</svg>")]

ROWS = ["Promoterless-donor-plasmid.svg",
        "Genomic-landing-pad---before-recombination.svg",
        "Integrated-cassette-in-genome---hybrid-sites.svg"]
SMALL = ["Paired-curved-site-exchange-strokes.svg",
         "Paired-curved-site-exchange-strokes-2.svg",
         "Integration-direction.svg", "Integration-arrowhead.svg",
         "Bxb1-recombination-label.svg"]
parts = {n: part(n) for n in ROWS + SMALL}

# where the donor's two att sites sit, so the crosses line up with them
donor = (d / ROWS[0]).read_text()
site_x = []
for m in re.finditer(r'id="Single_site_exterior[^"]*"[^>]*\sd="([^"]+)"', donor):
    xs = [float(v) for v in re.findall(r'[ML]\s*(-?[\d.]+)', m.group(1))]
    site_x.append((min(xs) + max(xs)) / 2)
site_x.sort()

K, PAD, GAP = 4.4, 16, 26
maxw = max(parts[n][0][2] for n in ROWS)
W = maxw * K + 2 * PAD
H = sum(parts[n][0][3] for n in ROWS) * K + 2 * PAD + GAP * 2 + 40
o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f">' % (W, H, W, H),
     '<rect width="%.0f" height="%.0f" fill="#fff"/>' % (W, H)]
for i, n in enumerate(ROWS + SMALL):
    o.append("<style>" + parts[n][1].replace(".cls-", ".p%d-cls-" % i) + "</style>")

def place(idx, name, x, y, k=K):
    vb, _, body = parts[name]
    b = body.replace('class="cls-', 'class="p%d-cls-' % idx)
    return ('<g transform="translate(%.1f,%.1f) scale(%s) translate(%s,%s)">%s</g>'
            % (x, y, k, -vb[0], -vb[1], b))

y = PAD
row_y = {}
for i, n in enumerate(ROWS):
    vb = parts[n][0]
    x = PAD + (maxw - vb[2]) * K / 2
    row_y[n] = (x, y)
    o.append(place(i, n, x, y))
    y += vb[3] * K + GAP

# the crosses, between the donor and the landing pad, on each att site
xd, yd = row_y[ROWS[0]]
cross_y = yd + parts[ROWS[0]][0][3] * K + 2
for j, sx in enumerate(site_x):
    nm = SMALL[j]
    vb = parts[nm][0]
    o.append(place(len(ROWS) + j, nm, xd + sx * K - vb[2] * K / 2, cross_y))

# The +Bxb1 arrow between the landing pad and the integrated locus. The two
# pieces are a bare shaft and a bare head, so the head has to be put at the
# shaft's far end: its apex is the bottom of its own box, and the shaft's line
# does not fill its box (0.2 of padding each way), so both are measured rather
# than assumed. Placing the head a fixed fraction down the shaft is what left
# the line poking out through it.
#
# The head also reaches HEAD_DROP past where the shaft ended, and the shaft is
# trimmed back to the head's base. Ending the shaft at the apex is not enough:
# its cap is round, so half a stroke still shows at the one point where the
# triangle has no width left to cover it.
SHAFT, HEAD = "Integration-direction.svg", "Integration-arrowhead.svg"
HEAD_DROP = 0.8
xi, yi = row_y[ROWS[2]]
ay = yi - GAP - 6
shaft_end = max(float(v) for v in
                re.findall(r'[ML]\s*-?[\d.]+\s+(-?[\d.]+)', parts[SHAFT][2]))
head_top = shaft_end + HEAD_DROP - parts[HEAD][0][3]
vb_s, st_s, body_s = parts[SHAFT]
parts[SHAFT] = (vb_s, st_s,
                re.sub(r'(\sd="M[\d.-]+\s+[\d.-]+\s+L[\d.-]+\s+)[\d.-]+',
                       lambda m: m.group(1) + "%g" % head_top, body_s))
for j, nm in enumerate([SHAFT, HEAD]):
    vb = parts[nm][0]
    o.append(place(len(ROWS) + 2 + j, nm, PAD + maxw * K / 2 - vb[2] * K / 2,
                   ay + (0 if nm == SHAFT else head_top * K)))
# the label sits beside the head, its baseline on the apex
vb = parts["Bxb1-recombination-label.svg"][0]
o.append(place(len(ROWS) + 4, "Bxb1-recombination-label.svg",
               PAD + maxw * K / 2 + 10,
               ay + (shaft_end + HEAD_DROP - 5.5) * K))
o.append("</svg>")
(d / "cassette_integration.svg").write_text("".join(o))
print("%.0fx%.0f  sites at %s" % (W, H, [round(v, 1) for v in site_x]))
