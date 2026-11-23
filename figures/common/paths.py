"""Central path handling. Import from every figure script."""
from pathlib import Path
import yaml

# Archive note: in the figures repository this tree sat one level deeper, so ROOT
# was parents[2]. Here it is figures/, hence parents[1]. This is the only edit
# made to the figure code when it was copied into the archive.
ROOT = Path(__file__).resolve().parents[1]
_cfg = yaml.safe_load((ROOT / "config" / "paths.yaml").read_text())

INPUTS = ROOT / _cfg["inputs"]
OUTPUT = ROOT / _cfg["output"]
FINAL = ROOT / _cfg["final"]
EULER_ROOT = _cfg["euler_root"]


def fig_output(name: str) -> Path:
    """Output folder for one figure, created if missing."""
    p = OUTPUT / name
    p.mkdir(parents=True, exist_ok=True)
    return p
