# rna_binder campaign control GUI

Local Streamlit app on your Mac that reads campaign state via eulermount
and runs commands on Euler over SSH.

## Setup (one-time)

```bash
# 1. Have eulermount mounted: ~/eulermount/lvalsecchi/rna_binder should exist
# 2. Make sure `ssh euler` works without a password (you already have this)
# 3. Install streamlit
cd ~/rna_binder_gui  # wherever you put these files
uv venv
source .venv/bin/activate
uv pip install streamlit pandas
```

If you don't want `uv`:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install streamlit pandas
```

## Running

```bash
source .venv/bin/activate
streamlit run app.py
```

Opens http://localhost:8501 in your browser.

To open in a chrome-less window instead of a tab:
```bash
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
  --app=http://localhost:8501
```

## Files

- `app.py` — Streamlit entry point
- `euler.py` — SSH wrappers (squeue, sbatch, scancel)
- `inspect_runs.py` — filesystem walker via eulermount
- `stacks.py` — per-stack metadata: which params can be overridden, defaults

## v1 (this version)

Inspect-only. Shows queue + recent runs + stack configs. Refresh button
clears the cache and re-fetches. No launch/cancel yet — those come in v2.
