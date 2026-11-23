# Figure and supplementary-table code

The code that draws every figure in the thesis and generates every supplementary
table. It was developed in a separate repository and copied here at the freeze so
that the deposit covers the presentation of the results as well as the pipeline
that produced them.

## Layout

```
Fig<N><letter>_<slug>/   one panel: make_figure.py, often a README.md explaining
                         the design decisions, sometimes an export_data.py that
                         writes the panel's source table
Fig<N>_assembly/         composes that figure's panels onto the A4 page
supp/<slug>/             the supplementary figures, same convention
supp/tables/             generators for Supplementary Tables (T3-T6 are generated
                         from primary records; T2 is a transcription rendered by
                         make_T2_rendering.py. Table S1 is not in this archive)
common/                  shared style: palette, glyphs, cartoon helpers, page
                         assembly, the Csy4 geometry, and paths.py
cartoon_library/         the structural cartoons
config/paths.yaml        where inputs and outputs live
environment/             requirements
run_all.sh               rebuild every panel
```

49 figure directories, 21 supplementary-figure directories, 191 Python files.

## Rebuilding

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r figures/environment/requirements.txt
PY=.venv/bin/python3 bash figures/run_all.sh
```

`run_all.sh` runs every `make_*.py` under `Fig*/` and `supp/*/`, skipping templates.
Individual panels run on their own: `python3 figures/Fig4c_fc_vs_splicing/make_figure.py`.
`common/STYLE_NOTES.md` is the style contract every panel follows, and
`common/CHIMERAX_RENDERING.md` covers the structure renders.

## What this tree does not contain

- **The rendered figures.** 318 MB of PDF and SVG, and every one of them is
  reproducible from this code. Scripts write into `output/`, which is created on
  first run.
- **The data the panels read.** Campaign tables, the design and re-prediction
  structures, and the screen read-outs are not in this repository; see
  `../provenance/data_availability.md` for where each tier lives. Several scripts
  also read directly from the cluster path in `config/paths.yaml`
  (`euler_root`), which only resolves inside ETH Zurich.
- **Working material**, left out deliberately when the tree was copied: the
  588 KB development log, two table-provenance notes, a superseded draft of the
  construct inventory, and five run logs.

## The one modification

`common/paths.py` resolved the repository root as `parents[2]`, which was correct
when this tree sat one level deeper. Here it is `parents[1]`. That single line is
the only change made to the figure code on the way in; everything else is byte for
byte what produced the figures.
