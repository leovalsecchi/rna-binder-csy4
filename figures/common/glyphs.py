# -*- coding: utf-8 -*-
import math

def _rot(pts, ang, cx, cy):
    a=math.radians(ang); c,s=math.cos(a), math.sin(a)
    return [((x-cx)*c-(y-cy)*s+cx, (x-cx)*s+(y-cy)*c+cy) for x,y in pts]

def _faces(cx, cy, turns, R, pitch, W, ang, n=28, chunk=math.pi):
    k = pitch/(2*math.pi); b = math.hypot(R,k); c=(W/2.0)/b
    T = 2*math.pi*turns; total = pitch*turns
    def edges(tt):
        x=k*tt
        return ((x-c*R, R*math.cos(tt)-c*k*math.sin(tt)),
                (x+c*R, R*math.cos(tt)+c*k*math.sin(tt)))
    out=[]; i=0
    while i*chunk < T-1e-9:
        t0=i*chunk; t1=min((i+1)*chunk, T)
        ts=[t0+(t1-t0)*j/n for j in range(n+1)]
        up=[edges(x)[0] for x in ts]; dn=[edges(x)[1] for x in ts]
        pts=up+list(reversed(dn))
        pts=[(x-total/2.0, y) for x,y in pts]
        pts=_rot(pts, ang, 0, 0)
        pts=[(x+cx, y+cy) for x,y in pts]
        out.append((math.sin((t0+t1)/2.0)>0, pts))
        i+=1
    return out

def _poly(pts, op):
    d="M"+" L".join(f"{x:.2f},{y:.2f}" for x,y in pts)+" Z"
    return f'<path d="{d}" fill="currentColor" fill-opacity="{op}" stroke="none"/>'

def ribbon_helix(cx=24, cy=24, turns=2.5, R=6.2, pitch=15.0, W=7.4, ang=-38,
                 front_op=0.62, back_op=1.0, n=28):
    """Two-tone ribbon alpha-helix from real helix geometry."""
    f=_faces(cx,cy,turns,R,pitch,W,ang,n)
    back="".join(_poly(p, back_op) for isf,p in f if not isf)
    front="".join(_poly(p, front_op) for isf,p in f if isf)
    return back+front

def ribbon_dissolve(cx=24, cy=24, turns=4.0, R=5.4, pitch=12.0, W=6.6, ang=-38,
                    solid_from=2, front_op=0.62, back_op=1.0, seed=3,
                    dot_r=1.5, jitter=2.4):
    """Same helix, tail rendered as scattered dots of its own silhouette."""
    import random
    rnd=random.Random(seed)
    f=_faces(cx,cy,turns,R,pitch,W,ang,n=20)
    nf=len(f); out=[]; dots=[]
    for idx,(isf,pts) in enumerate(f):
        if idx>=solid_from:
            out.append((isf,pts)); continue
        half=len(pts)//2
        up=pts[:half]; dn=list(reversed(pts[half:]))
        frac=(solid_from-idx)/max(solid_from,1)
        amp=jitter*(0.45+0.95*frac)
        for j in range(0,half,4):
            for w in (0.16, 0.84):
                x=up[j][0]*(1-w)+dn[j][0]*w
                y=up[j][1]*(1-w)+dn[j][1]*w
                dots.append((x+rnd.uniform(-amp,amp), y+rnd.uniform(-amp,amp)))
    s="".join(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{dot_r}" fill="currentColor" stroke="none"/>'
              for x,y in dots)
    s+="".join(_poly(p,back_op) for isf,p in out if not isf)
    s+="".join(_poly(p,front_op) for isf,p in out if isf)
    return s

def twisted_hairpin(cx=24, col_op=0.45, cross_at=0.72, loop_h=10.0,
                    ybot=43.0, ytop=22.0, amp=7.2, loop_taper=0.32):
    """RNA stem-loop with a single crossover in the stem.

    `cross_at` is where along the stem the strands swap, 0 at the foot and 1 at
    the loop. Pushing it up past 0.5 keeps the rails near-parallel for a longer
    run and tightens the loop. The angle is warped as theta = pi * u**p with p
    chosen so theta hits pi/2 exactly at cross_at.
    """
    p_warp = math.log(0.5) / math.log(cross_at)

    def rail(sign, u):
        y = ybot - (ybot - ytop) * u
        # taper the stem width only near the top, so the loop closes tighter
        a = amp * (1.0 - loop_taper * u ** 3)
        return (cx + sign * a * math.cos(math.pi * (u ** p_warp)), y)

    def path(sign, u0, u1, steps=22):
        pts = [rail(sign, u0 + (u1 - u0) * i / steps) for i in range(steps + 1)]
        return "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts)

    s = []
    s.append(f'<path d="{path(1, 0.0, cross_at + 0.005)}"/>')      # strand A, behind
    s.append(f'<path d="{path(-1, 0.0, 1.0)}"/>')                 # strand B, over
    s.append(f'<path d="{path(1, cross_at - 0.005, 1.0)}"/>')      # strand A, in front
    for u in (0.10, 0.30, 0.50, 0.66):
        (x1, y1), (x2, y2) = rail(1, u), rail(-1, u)
        s.append(f'<path d="M{x1:.2f},{y1:.2f} L{x2:.2f},{y2:.2f}" '
                 f'stroke-opacity="{col_op}"/>')
    lx, ly = rail(1, 1.0)
    rx, ry = rail(-1, 1.0)
    # closing pair: the loop's nucleotides touching across the top of the stem
    s.append(f'<path d="M{min(lx,rx)+1.2:.2f},{ly:.2f} '
             f'L{max(lx,rx)-1.2:.2f},{ry:.2f}" stroke-opacity="{col_op}"/>')
    s.append(f'<path d="M{lx:.2f},{ly:.2f} C{lx:.2f},{ly-loop_h:.1f} '
             f'{rx:.2f},{ry-loop_h:.1f} {rx:.2f},{ry:.2f}"/>')
    return "".join(s)



def _separate(pts, min_d):
    """drop points that sit closer than min_d, so dots never fuse into blobs"""
    keep=[]
    for p in pts:
        if all((p[0]-q[0])**2+(p[1]-q[1])**2 >= min_d*min_d for q in keep):
            keep.append(p)
    return keep

def ribbon_cloud(cx=24, cy=24, turns=2.5, R=6.2, pitch=15.0, W=7.4, ang=-38,
                 seed=5, dot_r=1.45, jitter=0.6, step=6, edges_only=True):
    """The same helix as a sparse point cloud tracing its own silhouette."""
    import random
    rnd=random.Random(seed)
    f=_faces(cx,cy,turns,R,pitch,W,ang,n=12)
    ws = (0.06, 0.94) if edges_only else (0.06, 0.5, 0.94)
    pts=[]
    for isf,poly in f:
        half=len(poly)//2
        up=poly[:half]; dn=list(reversed(poly[half:]))
        for j in range(0,half,step):
            for w in ws:
                x=up[j][0]*(1-w)+dn[j][0]*w
                y=up[j][1]*(1-w)+dn[j][1]*w
                pts.append((x+rnd.uniform(-jitter,jitter), y+rnd.uniform(-jitter,jitter)))
    pts=_separate(pts, dot_r*1.55)
    return "".join(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{dot_r}" fill="currentColor" stroke="none"/>'
                   for x,y in pts)

def ribbon_fade(cx=24, cy=24, turns=3.0, R=5.4, pitch=12.0, W=6.6, ang=-38,
                front_op=0.62, back_op=1.0, fade_frac=0.25, seed=4,
                dot_r=1.45, jitter=2.8, overlap=0, step=4):
    """One helix: tail is a sparse cloud of its own silhouette, head is solid ribbon.
    The cloud bleeds `overlap` faces into the solid part so the join reads as a fade."""
    import random
    rnd=random.Random(seed)
    f=_faces(cx,cy,turns,R,pitch,W,ang,n=12)
    nf=len(f); cut=max(1,int(round(nf*fade_frac)))
    body=[]; dots=[]
    for idx,(isf,poly) in enumerate(f):
        if idx>=cut:
            body.append((isf,poly, front_op if isf else back_op))
        if idx < cut+overlap:
            half=len(poly)//2
            up=poly[:half]; dn=list(reversed(poly[half:]))
            a = min(1.0, idx/float(cut))      # 0 at tail, 1 at the join
            amp=jitter*((1.0-a)**1.6)+0.25
            op = 1.0 if idx < cut else 0.55
            for j in range(0,half,step):
                for w in (0.06, 0.94):
                    x=up[j][0]*(1-w)+dn[j][0]*w
                    y=up[j][1]*(1-w)+dn[j][1]*w
                    dots.append((x+rnd.uniform(-amp,amp), y+rnd.uniform(-amp,amp), op))
    s="".join(_poly_op(p,o) for isf,p,o in body if not isf)
    s+="".join(_poly_op(p,o) for isf,p,o in body if isf)
    dots=_separate(dots, dot_r*1.55)
    s+="".join(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{dot_r}" fill="currentColor" '
               f'fill-opacity="{o:.2f}" stroke="none"/>' for x,y,o in dots)
    return s

def _poly_op(pts, op):
    d="M"+" L".join(f"{x:.2f},{y:.2f}" for x,y in pts)+" Z"
    return f'<path d="{d}" fill="currentColor" fill-opacity="{op:.3f}" stroke="none"/>'


def ribbon_ghost(cx=24, cy=24, turns=3.0, R=5.2, pitch=11.0, W=6.4, ang=-38,
                 front_op=0.62, back_op=1.0, ghost=0.30, seed=9,
                 dot_r=1.5, halo_min=2.6, halo_max=7.4, step=6, n_face=12, sep=2.9):
    """Faint ribbon with a sparse cloud of points scattered all around it."""
    import random, math as _m
    rnd=random.Random(seed)
    f=_faces(cx,cy,turns,R,pitch,W,ang,n=n_face)
    s ="".join(_poly_op(p, back_op*ghost)  for isf,p in f if not isf)
    s+="".join(_poly_op(p, front_op*ghost) for isf,p in f if isf)
    pts=[]
    for isf,poly in f:
        half=len(poly)//2
        up=poly[:half]; dn=list(reversed(poly[half:]))
        for j in range(0,half,step):
            for w in (0.0, 1.0):
                bx=up[j][0]*(1-w)+dn[j][0]*w
                by=up[j][1]*(1-w)+dn[j][1]*w
                a=rnd.uniform(0, 2*_m.pi)
                d=rnd.uniform(halo_min, halo_max)
                pts.append((bx+d*_m.cos(a), by+d*_m.sin(a)))
    rnd.shuffle(pts)
    pts=_separate(pts, dot_r*sep)
    s+="".join(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{dot_r}" fill="currentColor" stroke="none"/>'
               for x,y in pts)
    return s


def recycle_ring(cx=24, cy=24, r=16.5, a_start=104, sweep=312, sw=2.4, head=4.6, n=90):
    """Open circular arrow concentric with (cx,cy). Polyline so geometry is exact."""
    pts=[]
    for i in range(n+1):
        a=math.radians(a_start - sweep*i/n)
        pts.append((cx+r*math.cos(a), cy-r*math.sin(a)))
    d="M"+" L".join(f"{x:.2f},{y:.2f}" for x,y in pts)
    body=f'<path d="{d}" fill="none" stroke-width="{sw}"/>'
    ax,ay=pts[-1]
    bx,by=pts[-4]
    vx,vy=ax-bx, ay-by
    L=math.hypot(vx,vy) or 1.0
    vx,vy=vx/L, vy/L
    px,py=-vy,vx
    tipx,tipy = ax+vx*head*1.15, ay+vy*head*1.15
    p1=(ax+px*head*0.82, ay+py*head*0.82)
    p2=(ax-px*head*0.82, ay-py*head*0.82)
    body+=(f'<path d="M{tipx:.2f},{tipy:.2f} L{p1[0]:.2f},{p1[1]:.2f} '
           f'L{p2[0]:.2f},{p2[1]:.2f} Z" fill="currentColor" stroke="none"/>')
    return body


# ---------------------------------------------------------------- icon set --
# One style spec for every glyph, enforced on output:
#   viewBox 0 0 48 48 | stroke-width 2.4 | round caps and joins
#   stroke="currentColor" so colour comes from the parent | no transforms
#   content inside a 4 px margin | filled dots r=1.8
# Ribbon glyphs are fills (two tones from one colour via fill-opacity) and need
# roughly 80 px to read; the stroke glyphs survive down to 24 px.

SVG_HEAD = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 48 48" '
            'width="48" height="48" fill="none" stroke="currentColor" '
            'stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round">')
SVG_FOOT = '</svg>'


def dot(x, y, r=1.8, op=1.0):
    o = '' if op >= 0.999 else f' fill-opacity="{op}"'
    return f'<circle cx="{x}" cy="{y}" r="{r}" fill="currentColor" stroke="none"{o}/>'


# Residue identity is shown as a brightness pattern rather than a hue set, so the
# glyph still takes its colour from the parent. The same pattern is used in
# inverse-folding and refolding: one sequence written, then read back.
SEQ_X = (10, 17, 24, 31, 38)
SEQ_OP = (1.0, 0.55, 0.80, 0.95, 0.65)


def sequence_row(y, x0=10, x1=38):
    """Backbone line with five residues of varying brightness."""
    b = [f'<path d="M{x0},{y} L{x1},{y}" stroke-opacity="0.45"/>']
    b += [dot(x, y, op=o) for x, o in zip(SEQ_X, SEQ_OP)]
    return "".join(b)


def icon_tube():
    """Microcentrifuge tube: flat cap, straight walls, conical tip."""
    b = [
        # cap, a little wider than the body
        '<path d="M15.6,7.0 L32.4,7.0 A1.6,1.6 0 0 1 32.4,11.8 '
        'L15.6,11.8 A1.6,1.6 0 0 1 15.6,7.0 Z"/>',
        # body: parallel walls to the shoulder, then a conical tip
        '<path d="M17.6,11.8 L17.6,28.2 L23.0,40.4 '
        'A1.7,1.7 0 0 0 25.0,40.4 L30.4,28.2 L30.4,11.8"/>',
    ]
    # graduations, inner wall only
    for y in (16.6, 20.6, 24.6):
        b.append(f'<path d="M19.6,{y} L26.0,{y}" stroke-opacity="0.32"/>')
    # a little content near the tip
    for x, y, r, o in ((24.0, 33.2, 1.25, 0.70), (22.9, 36.3, 1.05, 0.50),
                       (25.2, 36.6, 1.05, 0.58)):
        b.append(dot(x, y, r, o))
    return "".join(b)


def icon_oligo_pool():
    """Oligo pool: short single strands of staggered length, one per binder.

    Brightness varies per strand (the SEQ_OP pattern) so the pool reads as
    many distinct sequences while the glyph keeps one parent colour.
    """
    b = []
    rows = ((9.0, 6.0, 36.0), (16.5, 11.0, 42.0), (24.0, 6.0, 38.0),
            (31.5, 10.0, 42.0), (39.0, 6.0, 34.0))
    for (y, x0, x1), op in zip(rows, SEQ_OP):
        n = 24
        pts = [(x0 + (x1 - x0) * i / n, y + 2.0 * math.sin((x0 + (x1 - x0) * i / n) * 2 * math.pi / 10.0))
               for i in range(n + 1)]
        d = "M" + " L".join(f"{x:.2f},{yy:.2f}" for x, yy in pts)
        b.append(f'<path d="{d}" stroke-opacity="{op}"/>')
    return "".join(b)


def _arrow_dn(x, y1, y2, head=3.3):
    """Short downward arrow between two glyph rows."""
    return (f'<path d="M{x},{y1} L{x},{y2}"/>'
            f'<path d="M{x-head},{y2-head} L{x},{y2+0.5} L{x+head},{y2-head}"/>')


def icon_rna_hairpin():
    """RNA stem-loop with a single crossover in the stem."""
    return twisted_hairpin()


def icon_protein():
    """Alpha-helix as a two-tone ribbon."""
    return ribbon_helix(24, 24, turns=2.5, R=6.2, pitch=15.0, W=7.4)


def icon_diffusion():
    """Faint ribbon inside a sparse surrounding point cloud."""
    return ribbon_ghost(24, 24, turns=3.0, R=5.0, pitch=10.4, W=6.2,
                        ghost=0.36, step=6)


def icon_hallucination():
    """Iterative optimisation: ring concentric with a centred ribbon."""
    return (recycle_ring(24, 24, r=16.5)
            + ribbon_helix(24, 24, turns=2.0, R=4.2, pitch=9.6, W=5.4))


def icon_inverse_folding():
    """Structure down to sequence."""
    return (ribbon_helix(24, 14, turns=2.0, R=4.8, pitch=10.5, W=6.0)
            + _arrow_dn(24, 27, 33)
            + sequence_row(41))


def icon_refolding():
    """Sequence back up to structure."""
    return (sequence_row(6)
            + _arrow_dn(24, 11, 17)
            + ribbon_helix(24, 34, turns=2.0, R=4.8, pitch=10.5, W=6.0))



def icon_binder_rna():
    """Binder design: a designed helix packed against the RNA hairpin."""
    return (twisted_hairpin(cx=15, ybot=42.0, ytop=19.0, amp=6.2, loop_h=8.5,
                            loop_taper=0.28)
            + ribbon_helix(33, 26, turns=2.0, R=5.0, pitch=10.5, W=6.0))


def icon_funnel():
    """In-silico filtering: many candidates in, few out."""
    b = ['<path d="M8.5,11 L39.5,11 L27,25.5 L27,35 L21,31.5 L21,25.5 Z"/>']
    b += [dot(x, y, 1.6) for x, y in ((14.5, 6.5), (24, 5.6), (33.5, 6.5))]
    b.append(dot(24, 41.5, 2.0))
    return "".join(b)


def icon_integration():
    """Site-specific integration: a donor cassette landing in the genome."""
    return (
        '<rect x="15" y="6" width="18" height="8.5" rx="2"/>'      # donor cassette
        '<path d="M24,16.5 L24,22.5"/>'
        '<path d="M21,19.8 L24,23.4 L27,19.8"/>'                   # arrow into locus
        '<path d="M5,31 L17,31"/><path d="M31,31 L43,31"/>'        # genome arms
        '<rect x="17" y="26.5" width="14" height="9" rx="2"/>'     # integrated
    )


def icon_long_reads():
    """Long-read sequencing: reads spanning the construct, staggered as aligned."""
    return (
        '<path d="M6,11 L42,11"/>'                                  # the construct
        '<path d="M6,7.5 L6,14.5"/><path d="M42,7.5 L42,14.5"/>'    # its two ends
        '<path d="M7,22 L31,22"/>'                                  # reads
        '<path d="M15,30 L41,30"/>'
        '<path d="M5,38 L25,38"/>'
    )


def icon_cytometry():
    """Flow cytometry: cells in single file past the interrogation point."""
    b = [
        '<path d="M13,5 C13,17 18,19 18,43"/>',          # hydrodynamic focusing
        '<path d="M35,5 C35,17 30,19 30,43"/>',
        '<path d="M4,23 L13,23"/>',                       # excitation
        '<path d="M9.8,20.4 L13.6,23 L9.8,25.6"/>',
        '<path d="M35,23 L44,23"/>',                      # emission to detector
        '<path d="M40.2,20.4 L44,23 L40.2,25.6"/>',
    ]
    b += [f'<circle cx="24" cy="{y}" r="2.8"/>' for y in (12, 23, 34)]
    return "".join(b)


def _wavy(y, x0, x1, amp=2.0, lam=10.0, n=24, op=None):
    """A single-stranded wave, the shape oligo-pool uses for RNA/ssDNA."""
    pts = [(x0 + (x1 - x0) * i / n,
            y + amp * math.sin((x0 + (x1 - x0) * i / n) * 2 * math.pi / lam))
           for i in range(n + 1)]
    d = "M" + " L".join(f"{x:.2f},{yy:.2f}" for x, yy in pts)
    o = "" if op is None else f' stroke-opacity="{op}"'
    return f'<path d="{d}"{o}/>'


def _duplex(y, x0, x1, gap=5.4, rungs=5):
    """Double-stranded DNA: two rails with base pairs between them."""
    b = [f'<path d="M{x0},{y - gap / 2:.2f} L{x1},{y - gap / 2:.2f}"/>',
         f'<path d="M{x0},{y + gap / 2:.2f} L{x1},{y + gap / 2:.2f}"/>']
    for i in range(rungs):
        x = x0 + (x1 - x0) * (i + 0.5) / rungs
        b.append(f'<path d="M{x:.2f},{y - gap / 2:.2f} L{x:.2f},'
                 f'{y + gap / 2:.2f}" stroke-opacity="0.38"/>')
    return "".join(b)


def icon_rt_cdna():
    """Reverse transcription: a single RNA strand becomes a cDNA duplex."""
    return (_wavy(9.0, 7.0, 41.0)
            + _arrow_dn(24, 17.0, 27.5)
            + _duplex(37.0, 7.0, 41.0, rungs=4))


def icon_pcr():
    """PCR: one duplex in, many out. Amplification is the count, not the size."""
    # the outputs carry no base pairs: at this size a rung between rails only
    # fills the gap that says "double-stranded"
    b = [_duplex(8.5, 13.0, 35.0, rungs=3), _arrow_dn(24, 15.0, 24.0)]
    # the gap WITHIN a duplex must be clearly smaller than the gap BETWEEN
    # them, or four rails read as four strands instead of two duplexes
    for y in (29.5, 41.5):
        b.append(_duplex(y, 6.0, 42.0, gap=4.0, rungs=0))
    return "".join(b)


def icon_barcode():
    """Barcoding: a tagged end on the construct, one bar pattern per condition.

    The bars are the barcode; the rule under them is the fragment they are
    ligated to, so the glyph says WHERE the tag sits and not just "barcode".
    """
    b = []
    for x, w, op in ((13.4, 2.6, 1.0), (17.8, 1.3, 0.55), (20.8, 2.0, 0.85),
                     (24.6, 1.3, 0.5), (27.6, 2.6, 1.0), (32.0, 1.4, 0.6)):
        b.append(f'<path d="M{x:.2f},12.0 L{x:.2f},29.0" '
                 f'stroke-width="{w:.2f}" stroke-opacity="{op}"/>')
    b.append('<path d="M9.0,34.5 L39.0,34.5"/>')
    b.append('<path d="M31.5,34.5 L39.0,34.5" stroke-opacity="0.38"/>')
    b.append('<path d="M9.0,31.0 L9.0,38.0"/>')
    b.append('<path d="M39.0,31.0 L39.0,38.0"/>')
    return "".join(b)


def icon_nanopore():
    """ONT: a strand threading a pore in the membrane, read as it passes."""
    b = ['<path d="M4.0,21.0 L18.0,21.0"/>', '<path d="M30.0,21.0 L44.0,21.0"/>',
         '<path d="M4.0,27.0 L18.0,27.0"/>', '<path d="M30.0,27.0 L44.0,27.0"/>']
    # the pore: two walls through the bilayer
    b.append('<path d="M18.0,18.0 L18.0,30.0"/>')
    b.append('<path d="M30.0,18.0 L30.0,30.0"/>')
    # the strand, threading top to bottom through it
    b.append(_wavy(9.0, 10.0, 38.0, amp=2.4, lam=13.0))
    b.append('<path d="M24.0,11.0 L24.0,40.0" stroke-opacity="0.55"/>')
    b.append(dot(24.0, 24.0, 1.5, 0.9))
    return "".join(b)


def icon_cell_lysis():
    """Lysis: the cell is open on one side and its RNA is coming out.

    The membrane is one arc with a gap rather than a closed outline, because
    the gap IS the event - a closed cell with strands beside it would read as
    a cell next to some RNA.
    """
    b = ['<path d="M31.0,11.5 A14.5,14.5 0 1 0 31.0,36.5"/>']
    # the two torn ends, thrown clear of the opening
    b.append('<path d="M33.6,9.6 A5.0,5.0 0 0 1 37.6,13.0" stroke-opacity="0.55"/>')
    b.append('<path d="M37.6,35.0 A5.0,5.0 0 0 1 33.6,38.4" stroke-opacity="0.55"/>')
    # RNA spilling through the gap
    b.append(_wavy(18.5, 29.0, 45.0, amp=1.7, lam=8.0))
    b.append(_wavy(29.5, 29.0, 45.0, amp=1.7, lam=8.0, op=0.6))
    return "".join(b)


def icon_read_filter():
    """Minimum-read filtering: per-barcode counts against a cutoff.

    Counts as bars with a dashed threshold across them, the ones that fall
    short faded out. A funnel would do for "filtering" in the abstract, but
    this says WHAT is filtered and on what.
    """
    thr = 26.0
    b = [f'<path d="M5.0,{thr} L43.0,{thr}" stroke-opacity="0.5" '
         f'stroke-dasharray="3.2 2.6"/>']
    for x, top, op in ((9.5, 11.0, 1.0), (16.5, 19.0, 1.0), (23.5, 32.5, 0.28),
                       (30.5, 15.0, 1.0), (37.5, 35.5, 0.28)):
        b.append(f'<path d="M{x},42.0 L{x},{top}" stroke-opacity="{op}"/>')
    return "".join(b)


ICONS = {
    "rna-hairpin": icon_rna_hairpin,
    "protein": icon_protein,
    "diffusion": icon_diffusion,
    "hallucination": icon_hallucination,
    "inverse-folding": icon_inverse_folding,
    "refolding": icon_refolding,
    "tube": icon_tube,
    "oligo-pool": icon_oligo_pool,
    "binder-rna": icon_binder_rna,
    "funnel": icon_funnel,
    "integration": icon_integration,
    "long-reads": icon_long_reads,
    "cytometry": icon_cytometry,
    "rt-cdna": icon_rt_cdna,
    "pcr": icon_pcr,
    "barcode": icon_barcode,
    "nanopore": icon_nanopore,
    "cell-lysis": icon_cell_lysis,
    "read-filter": icon_read_filter,
}


def icon_svg(name):
    """Standalone 48x48 SVG for one icon."""
    return SVG_HEAD + ICONS[name]() + SVG_FOOT


def icon_group(name, cx, cy, colour, k=2.3, sw=2.4):
    """The same glyph placed in a larger figure, scaled and coloured."""
    h = 24 * k
    return (f'<g transform="translate({cx-h},{cy-h}) scale({k})" fill="none" '
            f'stroke="{colour}" color="{colour}" stroke-width="{sw}" '
            f'stroke-linecap="round" stroke-linejoin="round">{ICONS[name]()}</g>')


def write_icons(outdir):
    """Write every icon as its own SVG. Returns the paths written."""
    from pathlib import Path
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    # Drop icons from a previous run whose glyph no longer exists, so the folder
    # never offers a name the figures cannot use.
    for stale in outdir.glob("*.svg"):
        if stale.stem not in ICONS:
            stale.unlink()
    written = []
    for name in ICONS:
        p = outdir / f"{name}.svg"
        p.write_text(icon_svg(name), encoding="utf-8")
        written.append(p)
    return written
