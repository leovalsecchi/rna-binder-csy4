Create env: `python3 -m venv .venv && source .venv/bin/activate && pip install -r environment/requirements.txt`
Or with conda: `conda create -n thesis_figs python=3.11 --file environment/requirements.txt`
Record any extra tools (PyMOL, ChimeraX versions) here.

On Leo's Mac the working interpreter is `~/.venvs/binder_map/bin/python3`
(Python 3.10, matplotlib 3.10.9). Rebuild everything with:

    PY=~/.venvs/binder_map/bin/python3 bash code/run_all.sh
