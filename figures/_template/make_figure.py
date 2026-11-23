"""Template: copy this folder to code/figNN_<name>/ and edit."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
style.apply()
out = paths.fig_output(FIG)

# --- load data from paths.INPUTS here ---

fig, ax = plt.subplots(figsize=(style.mm(style.CFG["width_single_mm"]), style.mm(60)))
ax.plot([0, 1], [0, 1], color=palette.MAIN)                 # single series: theme main
# per stack: ax.plot(x, y, color=palette.stack_color("BG+BIF"), label="BoltzGen + BoltzIF")
ax.set_xlabel("x")
ax.set_ylabel("y")
style.save(fig, out, "panel_a")
print("wrote", out)
