"""
rna_binder campaign control GUI — v1 (inspect-only).

Usage:
    streamlit run app.py

Shows the SLURM queue + configured stacks, with launch + cancel controls.
"""
from __future__ import annotations
import os
import json
import streamlit as st
import pandas as pd
import yaml

import euler
import inspect_runs as inspect
import analysis_ui
from stacks import (
    STACKS, StackSpec, ParamSpec, REFOLD_AF3_PARAMS, REFOLD_B2_PARAMS, GPU_GPUMEM,
    REFOLD_ESM_PARAMS, REFOLD_ESM_CONFIG,
    REFOLD_AF3_CONFIG, REFOLD_B2_CONFIG, BULK_PARTITION, REFOLD_GPU_STRATEGIES,
    ALL_24G_NODES, AMPERE_PLUS_NODES, estimate_designs,
)


def pair_gpumem(overrides: list[tuple[str, str]]) -> list[tuple[str, str]]:
    """Pair a slurm.gpus override with a matching slurm.gres (gpumem) override.

    Euler's gpumem gres is a minimum-VRAM filter, so it must fit the chosen
    card; without this, picking a 24 GB card under a config that asks for
    gpumem:40g yields 'Requested node configuration is not available'.
    """
    out: list[tuple[str, str]] = []
    for path, val in overrides:
        if path == "slurm.gpus" and val == "turing_24g":
            # "All usable 24 GB" pool: 3090 + 4090 + titan + quadro (ALL_24G_NODES).
            # gpus=1 + the node allow-list, so SLURM grabs whichever 24 GB card is
            # free first (most convenient / shortest queue). Pin by NODELIST, not
            # GPU type — on Euler a typed --gpus=<turing>:1 routes loosely (a quadro
            # request landed on a titan node). Generators run on Turing sm_75 with
            # no flag (GENERATOR_GPU_COMPAT.md). Stays in gpuhe (no partition change).
            out.append(("slurm.gpus", "1"))
            out.append(("slurm.gres", "gpumem:20g"))
            out.append(("slurm.nodelist", ALL_24G_NODES))
        elif path == "slurm.gpus" and val == "ampere_24g":
            # Ampere/Ada 24 GB pool: 3090 + 4090 only (AMPERE_PLUS_NODES, no Turing).
            # Required for HalluDesign (Protenix Evoformer needs sm_80+; Turing sm_75 fails
            # with "Unsupported GPU and data type combination" -> RNA-only designs).
            out.append(("slurm.gpus", "1"))
            out.append(("slurm.gres", "gpumem:20g"))
            out.append(("slurm.nodelist", AMPERE_PLUS_NODES))
        elif path == "slurm.gpus" and val in GPU_GPUMEM:
            out.append((path, val))
            out.append(("slurm.gres", GPU_GPUMEM[val]))
        else:
            out.append((path, val))
    return out


# Map the refold walltime dropdown to the gpupr partition tier (a100 strategies).
_GPUPR_TIER = {"04:00:00": "4h", "24:00:00": "24h", "120:00:00": "120h"}


def expand_refold_gpu(overrides: list[tuple[str, str]]) -> list[tuple[str, str]]:
    """Expand a refold GPU *strategy* into its full SLURM override set.

    The refold GPU dropdown stores a strategy key in slurm.gpus (e.g.
    'ampere_24g'); turn it into slurm.gpus=1 + slurm.gres=gpumem:.. + optional
    slurm.nodelist (+ slurm.partition=gpupr.<tier> for the a100 strategies, since
    a100 lives in gpupr not gpuhe). Supervisor-style VRAM floor + node-range
    allow-list; mirrors refold__{af3,b2} default.yaml; see
    stacks.REFOLD_GPU_STRATEGIES. The per-type pair_gpumem() path stays for the
    generation stacks.
    """
    time_val = next((v for p, v in overrides if p == "slurm.time"), "04:00:00")
    out: list[tuple[str, str]] = []
    for path, val in overrides:
        if path == "slurm.gpus" and val in REFOLD_GPU_STRATEGIES:
            _label, gres, nodelist, partition = REFOLD_GPU_STRATEGIES[val]
            out.append(("slurm.gpus", "1"))
            out.append(("slurm.gres", gres))
            # Always emit nodelist so it overrides the YAML default: a value pins
            # the node range; "" clears it (e.g. the a100 strategies, which select
            # via the gpupr partition + gpumem floor and must NOT carry the
            # gpuhe-only eu-g4/eu-g6 default into gpupr).
            out.append(("slurm.nodelist", nodelist))
            if partition == "gpupr":
                out.append(("slurm.partition", f"gpupr.{_GPUPR_TIER.get(time_val, '4h')}"))
            elif partition == "gpu":
                # public gpu.* pool, tier picked by walltime (gpu.4h/24h/120h). The
                # gpumem:8g floor (from the strategy) keeps it off wicky gpuhe.
                out.append(("slurm.partition", f"gpu.{_GPUPR_TIER.get(time_val, '4h')}"))
        else:
            out.append((path, val))
    return out


def force_partition(overrides: list[tuple[str, str]], partition: str) -> list[tuple[str, str]]:
    """Override slurm.partition in an override list (for the 'run on bulk' toggles).

    Replaces an existing slurm.partition entry in place, else appends one. Leaves
    slurm.gpus/gres untouched — `-p bulk` maps to gpuhe.bulk, the same high-end
    GPUs, just preemptible. Safe because every generator + refold sets --requeue.
    """
    out, seen = [], False
    for path, val in overrides:
        if path == "slurm.partition":
            out.append((path, partition)); seen = True
        else:
            out.append((path, val))
    if not seen:
        out.append(("slurm.partition", partition))
    return out


# ─── shared helpers ──────────────────────────────────────────────────────────


def _coerce(x: str):
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return x


@st.cache_data(show_spinner=False)
def _load_config_yaml(abs_path: str) -> dict:
    """Load a stack/refold default.yaml off the mount so widgets can be seeded
    from the live cluster config. Cached (the sidebar Refresh clears it). Returns
    {} on any error, so seeding silently falls back to ParamSpec.default — the
    GUI still works when the eulermount is down."""
    try:
        with open(abs_path) as f:
            return yaml.safe_load(f) or {}
    except (OSError, yaml.YAMLError):
        return {}


def _yaml_value(cfg: dict, dotted: str):
    """Resolve a dotted path (e.g. 'generator_params.mpnn_cli.batch_size') in a
    loaded config dict, or None if any segment is missing."""
    cur = cfg
    for k in dotted.split("."):
        if isinstance(cur, dict) and k in cur:
            cur = cur[k]
        else:
            return None
    return cur


def _seed_for(p: ParamSpec, cfg: dict) -> str:
    """Widget seed for p: the live YAML value (as a widget-ready string), else
    p.default. Lists (e.g. model_seeds) become the comma string the as_list
    widget expects; YAML bools stringify to 'True'/'False' which render_param
    lowercases for the checkbox."""
    val = _yaml_value(cfg, p.path)
    if val is None:
        return p.default
    if p.as_list and isinstance(val, list):
        return ",".join(str(x) for x in val)
    return str(val)


def render_param(p: ParamSpec, key_prefix: str, seed: str | None = None) -> tuple[str, str] | None:
    """Render one override widget; return (dotted_path, value_string), or None for
    info-only widgets (rendered for display but contribute no override).

    seed, when given, is the live-YAML value used as the widget's initial value
    (falling back to p.default). It only affects the first render of each widget;
    user edits then persist via session_state."""
    key = f"{key_prefix}_{p.path}"
    seed_val = p.default if seed is None else seed
    # Info box = the help text plus, when present, the repo/paper-default citation.
    help_text = p.help
    if p.paper:
        prefix = f"{help_text}\n\n" if help_text else ""
        help_text = f"{prefix}📄 **Repo/paper default:** {p.paper}"
    if p.is_bool:
        checked = st.checkbox(p.name, value=str(seed_val).lower() == "true",
                              key=key, help=help_text or None, disabled=p.disabled)
        val = "true" if checked else "false"
    elif p.choices:
        # Seed from the live YAML when it's a valid choice; else fall back to the
        # ParamSpec default (e.g. a legacy slurm.time like 12:00:00 that isn't one
        # of the new tier choices), else the first option.
        if seed_val in p.choices:
            idx = p.choices.index(seed_val)
        elif p.default in p.choices:
            idx = p.choices.index(p.default)
        else:
            idx = 0
        fmt = (lambda v: p.choice_labels.get(v, v)) if p.choice_labels else (lambda v: v)
        val = st.selectbox(p.name, p.choices, index=idx, key=key,
                           help=help_text or None, format_func=fmt, disabled=p.disabled)
    else:
        val = st.text_input(p.name, value=seed_val, key=key, help=help_text or None,
                            disabled=p.disabled)
    if p.as_list:
        items = [x.strip() for x in str(val).split(",") if x.strip()]
        val = json.dumps([_coerce(x) for x in items])
    if p.info_only:
        return None
    return (p.path, val)


# ─── page setup ────────────────────────────────────────────────────────────
_FAVICON = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets", "favicon.png")
st.set_page_config(
    page_title="rna_binder",
    page_icon=_FAVICON if os.path.exists(_FAVICON) else "🧬",
    layout="wide",
)
st.markdown(
    "<style>"
    ".block-container{padding-top:1rem;padding-bottom:1rem;max-width:1200px;}"
    "</style>",
    unsafe_allow_html=True,
)

# ─── sidebar: config ───────────────────────────────────────────────────────
default_mount = os.path.expanduser(inspect.EULERMOUNT_ROOT_DEFAULT)
mount = st.sidebar.text_input("eulermount root", value=default_mount)
target = st.sidebar.text_input("Target", value="csy4_crrna_free")
user = st.sidebar.text_input("SLURM user", value="lvalsecchi")
remote_cwd = st.sidebar.text_input(
    "remote project dir",
    value="/cluster/work/wicky/lvalsecchi/rna_binder",
    help="Working dir on Euler where launch.py runs (launch.py is cd'd here).",
)
if st.sidebar.button("Refresh", use_container_width=True):
    st.cache_data.clear()
    st.rerun()


# ─── data fetchers (cached for the page lifetime; click Refresh to clear) ──
@st.cache_data(show_spinner=False)
def fetch_queue(user: str):
    try:
        return euler.squeue(user)
    except Exception as e:
        st.error(f"squeue failed: {e}")
        return []


# ─── main ──────────────────────────────────────────────────────────────────
st.title("rna_binder campaign control")
st.caption(f"target: `{target}`  •  mount: `{mount}`")

# Top-level view switch. Analysis renders its own page and stops, so the
# launch/monitor body below only runs for the default view (no re-indent needed).
view = st.radio("view", ["Launch & monitor", "Analysis"], horizontal=True,
                label_visibility="collapsed")
if view == "Analysis":
    analysis_ui.render(mount)
    st.stop()

st.subheader("SLURM queue")
jobs = fetch_queue(user)
if not jobs:
    st.info("No jobs running.")
else:
    qdf = pd.DataFrame([{
        "jobid": j.jobid, "name": j.name, "state": j.state,
        "elapsed": j.elapsed, "nodelist/reason": j.reason,
    } for j in jobs])
    st.dataframe(qdf, hide_index=True, use_container_width=True)

    # Cancel jobs straight from the queue (euler.scancel). Multiselect + button
    # rather than per-row buttons so a stray click can't kill a job.
    job_opts = {f"{j.jobid} — {j.name} [{j.state}]": j.jobid for j in jobs}
    qc1, qc2 = st.columns([4, 1])
    to_cancel = qc1.multiselect(
        "Cancel jobs", list(job_opts), label_visibility="collapsed",
        placeholder="Select job(s) to cancel…",
    )
    if qc2.button("Cancel", use_container_width=True, disabled=not to_cancel):
        for disp in to_cancel:
            try:
                euler.scancel(job_opts[disp])
            except Exception as e:  # ssh timeout / connection issues
                st.error(f"scancel {job_opts[disp]} failed: {e}")
        st.toast(f"Cancelled {len(to_cancel)} job(s)")
        st.cache_data.clear()  # re-fetch the queue without the killed jobs
        st.rerun()

st.subheader("Configured stacks")
sdf = pd.DataFrame([{
    "stack": s.key, "label": s.label, "config": s.config,
    "params": ", ".join(p.name for p in s.params),
} for s in STACKS])
st.dataframe(sdf, hide_index=True, use_container_width=True)


# ─── new campaign (v2: launch) ──────────────────────────────────────────────
st.divider()
st.subheader("New campaign")

label = st.text_input(
    "Run label", value="midrun_v2",
    help="Prefix for every run dir created by this launch.",
)


# ── presets ──────────────────────────────────────────────────────────────
# Two preset shapes coexist:
#   * Built-ins ({"num_runs", "params": {path: val}}) — path-based, applied to
#     every stack that exposes a matching ParamSpec.path. Easy to define;
#     stack-agnostic so e.g. slurm.time hits every stack at once.
#   * User-saved ({"num_runs", "stacks": {stack.key: {path: val}}}) — full
#     widget-state snapshots, restored exactly. Written by the Save preset
#     button to PRESETS_DIR.
# _apply_preset handles both shapes; snapshot_state writes the second.
PRESETS_DIR = "/Volumes/eulermount/lvalsecchi/rna_binder/_gui_presets"


def load_user_presets() -> dict[str, dict]:
    """Read every saved preset JSON from PRESETS_DIR.

    Returns {} silently if the eulermount is down — no toast spam on every
    page render. Bad JSON files are skipped, not surfaced.
    """
    out: dict[str, dict] = {}
    try:
        if not os.path.isdir(PRESETS_DIR):
            return out
        for fn in sorted(os.listdir(PRESETS_DIR)):
            # Skip macOS AppleDouble sidecars (._<name>) and other dotfiles;
            # AppleDouble files are binary and crash json.load with UnicodeDecodeError.
            if fn.startswith(".") or not fn.endswith(".json"):
                continue
            try:
                with open(os.path.join(PRESETS_DIR, fn)) as f:
                    out[fn[:-5]] = json.load(f)
            except (OSError, json.JSONDecodeError, UnicodeDecodeError):
                continue
    except OSError:
        pass
    return out


BUILTIN_PRESETS: dict[str, dict] = {
    "(none)": {},
    "Paper / repo defaults": {
        # Each method's published repo/paper default for the algorithm-QUALITY
        # knobs (see every field's ⓘ for the citation). SCALE knobs
        # (num_designs / n_batches / design_samples) and refold sample counts
        # are intentionally NOT set: those are budget/throughput choices, and
        # BoltzGen's documented default (10k backbones) would be a launch footgun.
        # num_runs is left untouched too.
        "params": {
            "generator_params.diffusion_batch_size": "8",      # RFD3 (foundry rfdiffusion3.yaml)
            "generator_params.params.num_seqs": "8",           # HalluDesign argparse default
            "generator_params.params.num_recycles": "15",      # HD examples 10–30; we use 15
            "generator_params.params.ref_time_steps": "60",    # HD live config (50 local / 150 design)
            "generator_params.params.mpnn_num_packs": "1",     # HD argparse default
            "generator_params.bif_cli.inverse_fold_num_sequences": "1",  # BoltzGen repo default
            "generator_params.mpnn_cli.batch_size": "8",       # RFdiffusion→MPNN convention: 8 seqs/backbone
            "generator_params.mpnn_cli.number_of_batches": "1",
            "generator_params.mpnn_cli.temperature": "0.1",    # ProteinMPNN/LigandMPNN standard
        },
    },
    "Smoke test": {
        "num_runs": 1,
        "params": {
            "generator_params.n_batches": "1",
            "generator_params.params.num_seqs": "2",
            "generator_params.params.num_recycles": "2",
            "generator_params.cli.flags.num_designs": "1",
            "generator_params.bg_cli.num_designs": "1",
            "generator_params.design_samples": "2",
            "slurm.time": "01:00:00",
        },
    },
    "Production (200 designs)": {
        "num_runs": 1,
        "params": {
            "generator_params.cli.flags.num_designs": "200",
            "generator_params.bg_cli.num_designs": "200",
            "generator_params.design_samples": "200",
            "generator_params.n_batches": "25",
            "generator_params.params.num_seqs": "10",
            "generator_params.params.num_recycles": "20",
            "generator_params.params.ref_time_steps": "30",
            "slurm.time": "23:30:00",
        },
    },
    "Production 1k (5×200)": {
        "num_runs": 5,
        "params": {
            "generator_params.cli.flags.num_designs": "200",
            "generator_params.bg_cli.num_designs": "200",
            "generator_params.design_samples": "200",
            "generator_params.n_batches": "25",
            "generator_params.params.num_seqs": "10",
            "generator_params.params.num_recycles": "20",
            "generator_params.params.ref_time_steps": "30",
            "slurm.time": "23:30:00",
        },
    },
}

# Merge user-saved presets at module-render time. Saved names win over
# built-ins on collision — gives the user an escape hatch to override a
# built-in. Order: built-ins first (familiar), then saved ones alphabetically.
PRESETS: dict[str, dict] = {**BUILTIN_PRESETS, **load_user_presets()}


def _apply_preset():
    """on_change callback for the preset dropdown.

    Writes preset values into the matching widget session_state keys BEFORE
    the widgets render this run, so the displayed values reflect the preset.
    Handles both built-in (path-based) and saved (stack-keyed) shapes.
    """
    p = PRESETS.get(st.session_state.get("preset_choice", "(none)"))
    if not p:
        return
    # Legacy/built-in single num_runs → seed every stack's per-stack widget.
    if "num_runs" in p:
        for stack in STACKS:
            st.session_state[f"num_runs_{stack.key}"] = int(p["num_runs"])
    # New per-stack shape (written by snapshot_state).
    if "num_runs_by_stack" in p:
        for stack_key, nr in p["num_runs_by_stack"].items():
            st.session_state[f"num_runs_{stack_key}"] = int(nr)
    # Built-in path-based: apply to every stack with a matching path.
    if "params" in p:
        params = p["params"]
        for stack in STACKS:
            for ps in stack.params:
                if ps.path in params:
                    st.session_state[f"param_{stack.key}_{ps.path}"] = str(params[ps.path])
    # User-saved stack-keyed: restore each widget exactly.
    if "stacks" in p:
        for stack_key, stack_vals in p["stacks"].items():
            for path, val in stack_vals.items():
                st.session_state[f"param_{stack_key}_{path}"] = str(val)


def snapshot_state() -> dict:
    """Capture the current form state as a stack-keyed preset.

    Reads st.session_state for each ParamSpec widget; falls back to the
    ParamSpec default for any widget the user hasn't touched. Used by the
    Save preset button.
    """
    return {
        "num_runs_by_stack": {
            stack.key: int(st.session_state.get(f"num_runs_{stack.key}", 1))
            for stack in STACKS
        },
        "stacks": {
            stack.key: {
                ps.path: st.session_state.get(f"param_{stack.key}_{ps.path}", ps.default)
                for ps in stack.params
            }
            for stack in STACKS
        },
    }


st.selectbox(
    "Load preset", list(PRESETS.keys()),
    key="preset_choice", on_change=_apply_preset,
    help="Bulk-set matching fields across every stack. Tweak individual "
         "values afterwards — presets just seed the form. 'Paper / repo "
         "defaults' fills each method's published values for the quality "
         "knobs only (scale & refold sample counts are left as-is).",
)

st.caption("**Num runs** is set per stack — each stack's card below has its own "
           "input, so you can choose independently. Each stack launches that "
           "many parallel gen jobs writing to a shared `<label>/designs/` dir.")

st.caption("Check the stacks to launch and adjust per-stack overrides. "
           "Unchecked stacks are skipped. Field defaults are seeded at page load "
           "from each stack's live default.yaml on the mount.")

N_COLS = 3
# Display labels for named ParamSpec groups (rendered inside collapsed expanders
# under the always-visible main params). Order here defines render order.
GROUP_LABELS = [("bif", "BoltzIF options"), ("limpnn", "LigandMPNN options"), ("slurm", "GPU / SLURM")]


def render_stack_params(stack: StackSpec, key_prefix: str) -> list[tuple[str, str]]:
    """Render a stack's ParamSpecs, grouping by ParamSpec.group.

    Default group ("") renders flat at the top; named groups render inside a
    collapsed st.expander each. Widget defaults are seeded from the stack's live
    default.yaml on the mount (falling back to ParamSpec.default per missing key),
    so the form reflects the current cluster config. Returns the override list in
    stacks.py order.
    """
    cfg = _load_config_yaml(os.path.join(mount, stack.config))
    by_group: dict[str, list[ParamSpec]] = {}
    for p in stack.params:
        by_group.setdefault(p.group or "", []).append(p)
    overrides: list[tuple[str, str]] = []
    bulk = False
    for p in by_group.get("", []):
        overrides.append(render_param(p, key_prefix, seed=_seed_for(p, cfg)))
    for group_key, group_label in GROUP_LABELS:
        if group_key in by_group:
            with st.expander(group_label, expanded=False):
                if group_key == "slurm":
                    st.caption(
                        "No partition picker: Euler auto-selects the gpuhe tier "
                        "from your slurm time. The only lever is the low-priority "
                        "**bulk** queue below."
                    )
                    bulk = st.checkbox(
                        "Run on bulk (low priority, preemptible)",
                        value=False, key=f"{key_prefix}_bulk",
                        help="Submit with `-p bulk` → gpuhe.bulk: runs on idle "
                             "high-end GPUs at low priority and can be preempted "
                             "(jobs set --requeue so they resume). Use for "
                             "non-urgent batches that shouldn't compete with your "
                             "other jobs.",
                    )
                for p in by_group[group_key]:
                    overrides.append(render_param(p, key_prefix, seed=_seed_for(p, cfg)))
    if bulk:
        overrides = force_partition(overrides, BULK_PARTITION)
    return overrides


# Collect (stack, overrides) for every checked stack. Overrides come straight
# from stacks.py paths — we never invent paths here.
selected: list[tuple[StackSpec, list[tuple[str, str]], int]] = []
for row_start in range(0, len(STACKS), N_COLS):
    cols = st.columns(N_COLS)
    for col, stack in zip(cols, STACKS[row_start:row_start + N_COLS]):
        with col, st.container(border=True):
            checked = st.checkbox(stack.label, value=False, key=f"chk_{stack.key}")
            stack_num_runs = st.number_input(
                "Num runs", min_value=1, value=1, step=1,
                key=f"num_runs_{stack.key}",
                help="Parallel independent gen jobs for THIS stack. When >1, all "
                     "runs write to a shared `<label>/designs/` dir.",
            )
            overrides = render_stack_params(stack, f"param_{stack.key}")
            if checked:
                selected.append((stack, pair_gpumem(overrides), int(stack_num_runs)))


def build_jobs(selected: list[tuple[StackSpec, list[tuple[str, str]]]]):
    """One submission per selected stack.

    Returns a list of (label, config, overrides) submission jobs. The BoltzGen
    multif stack already runs both inverse folders (BoltzIF + LigandMPNN) off a
    single shared backbone, so there is no longer any cross-stack merging here.
    """
    return [(stack.label, stack.config, overrides, num_runs)
            for stack, overrides, num_runs in selected]


auto_refold = st.checkbox(
    "Auto-refold", value=True,
    help="Queue AF3 + B2 refolds with an afterok dependency on each design job.",
)

# Global refold settings — applied to every auto-refold in this launch via
# launch.py's --refold-af3-override / --refold-b2-override flags.
refold_af3_overrides: list[tuple[str, str]] = []
refold_b2_overrides: list[tuple[str, str]] = []
refold_esm_overrides: list[tuple[str, str]] = []
refold_on_bulk = False
with st.expander("Refold settings (applied to all auto-refolds)",
                 expanded=False):
    if not auto_refold:
        st.caption("Enable **Auto-refold** above for these to take effect.")
    refold_on_bulk = st.checkbox(
        "Run on bulk (low priority, preemptible)",
        value=False, key="refold_on_bulk",
        help="Submit refolds with `-p bulk` → gpuhe.bulk (low priority, "
             "preemptible; --requeue is set so they resume).",
    )
    _af3_cfg = _load_config_yaml(os.path.join(mount, REFOLD_AF3_CONFIG))
    _b2_cfg = _load_config_yaml(os.path.join(mount, REFOLD_B2_CONFIG))
    _esm_cfg = _load_config_yaml(os.path.join(mount, REFOLD_ESM_CONFIG))
    ra, rb, rc = st.columns(3)
    with ra:
        st.markdown("**AlphaFold3**")
        for p in REFOLD_AF3_PARAMS:
            ov = render_param(p, "refold_af3", seed=_seed_for(p, _af3_cfg))
            if ov is not None:
                refold_af3_overrides.append(ov)
    with rb:
        st.markdown("**Boltz-2**")
        for p in REFOLD_B2_PARAMS:
            ov = render_param(p, "refold_b2", seed=_seed_for(p, _b2_cfg))
            if ov is not None:
                refold_b2_overrides.append(ov)
    with rc:
        st.markdown("**ESMFold2**")
        for p in REFOLD_ESM_PARAMS:
            ov = render_param(p, "refold_esm", seed=_seed_for(p, _esm_cfg))
            if ov is not None:
                refold_esm_overrides.append(ov)

refold_af3_overrides = expand_refold_gpu(refold_af3_overrides)
refold_b2_overrides = expand_refold_gpu(refold_b2_overrides)
refold_esm_overrides = expand_refold_gpu(refold_esm_overrides)
if refold_on_bulk:
    # a100 (gpupr) has no bulk tier, so leave a gpupr-pinned partition alone;
    # bulk only applies to the gpuhe-pool strategies.
    def _bulk_or_keep(ovs):
        if any(p == "slurm.partition" and str(v).startswith("gpupr") for p, v in ovs):
            return ovs
        return force_partition(ovs, BULK_PARTITION)
    refold_af3_overrides = _bulk_or_keep(refold_af3_overrides)
    refold_b2_overrides = _bulk_or_keep(refold_b2_overrides)
    refold_esm_overrides = _bulk_or_keep(refold_esm_overrides)

with st.expander("Preview launch command(s)", expanded=False):
    if not selected:
        st.caption("Select at least one stack to see the command(s).")
    else:
        for job_label, job_config, overrides, job_num_runs in build_jobs(selected):
            st.markdown(f"**{job_label}**  ·  {job_num_runs} run(s)")
            preview_cmd = euler.build_launch_cmd(
                remote_cwd, job_config, overrides, label, auto_refold,
                num_runs=int(job_num_runs),
                refold_af3_overrides=refold_af3_overrides,
                refold_b2_overrides=refold_b2_overrides,
                refold_esm_overrides=refold_esm_overrides,
            )
            st.code(preview_cmd, language="bash")

# Pre-launch scale readout: how many designs / jobs this launch will produce.
# Estimates only (see stacks.estimate_designs) — labelled as such.
if selected:
    lines, total, have_total, n_gen_jobs = [], 0, True, 0
    for stack, ov_list, nr in selected:
        n_gen_jobs += nr
        cnt, formula = estimate_designs(stack.key, dict(ov_list))
        if cnt is None:
            lines.append(f"- **{stack.label}**: (no estimate, × {nr} run(s))")
            have_total = False
            continue
        per = cnt * nr
        total += per
        run_suffix = f" × {nr} runs" if nr > 1 else ""
        lines.append(f"- **{stack.label}**: ≈{per:,} designs  ({formula}{run_suffix})")
    head = f"**Estimated scale** — {n_gen_jobs} generation job(s)"
    if auto_refold:
        head += " + AF3 & B2 refold per run"
    if have_total:
        head += f"  •  ≈{total:,} designs total"
    st.markdown(head + "\n" + "\n".join(lines))

preset_save_name = st.text_input(
    "Preset name (for Save preset)", value="", key="preset_save_name",
    help=f"Snapshot every current form value to {PRESETS_DIR}/<name>.json. "
         "Appears in 'Load preset' on the next page render. Overwrites if name exists.",
)

c1, c2, c3 = st.columns(3)
launch_clicked = c1.button("Launch selected", type="primary", use_container_width=True)
dryrun_clicked = c2.button("Dry-run", use_container_width=True)
save_clicked = c3.button("Save preset", use_container_width=True)

if save_clicked:
    name = (preset_save_name or "").strip()
    if not name:
        st.warning("Enter a preset name first.")
    elif "/" in name or "\\" in name or name.startswith("."):
        st.warning("Preset name cannot contain slashes or start with '.'")
    elif name in BUILTIN_PRESETS:
        st.warning(f"'{name}' is a built-in preset; pick a different name.")
    else:
        try:
            os.makedirs(PRESETS_DIR, exist_ok=True)
            path = os.path.join(PRESETS_DIR, f"{name}.json")
            with open(path, "w") as f:
                json.dump(snapshot_state(), f, indent=2)
            st.success(f"Saved preset '{name}' -> {path}")
            st.rerun()
        except OSError as e:
            st.error(f"Save failed (is eulermount up?): {e}")

# Results from a just-completed launch survive the auto-refresh rerun below.
if "launch_results" in st.session_state:
    for lbl, r in st.session_state.pop("launch_results"):
        if r["ok"] and r.get("detached"):
            st.toast(f"{lbl}: launch started in background")
            st.success(
                f"**{lbl}** — launch.py is running on Euler in the background. "
                "Jobs enter the SLURM queue as they submit (this can take a minute "
                "or two for a big multi-run launch) — hit **Refresh** to watch them. "
                f"\n\nLog: `{r.get('log_path') or '(see launch_logs/)'}`  "
                "— check it if jobs don't appear."
            )
        elif r["ok"]:
            st.toast(f"{lbl}: submitted jobid {r['jobid']}")
            extra = f"  (refolds: {', '.join(r['refolds'])})" if r["refolds"] else ""
            st.success(f"{lbl} -> jobid {r['jobid']}{extra}")
        else:
            st.error(f"{lbl} failed (exit {r['rc']})\n\n{r['err']}")

def _ready_to_submit() -> bool:
    """Common pre-flight: at least one stack must be selected."""
    if not selected:
        st.warning("No stacks selected.")
        return False
    return True


if dryrun_clicked and _ready_to_submit():
    for job_label, job_config, overrides, job_num_runs in build_jobs(selected):
        with st.expander(f"Dry-run: {job_label} ({job_num_runs} run(s))", expanded=True):
            with st.spinner("Running launch.py --dry-run on Euler..."):
                res = euler.sbatch_via_launcher(
                    remote_cwd, job_config, overrides, label, auto_refold, dry_run=True,
                    num_runs=int(job_num_runs),
                    refold_af3_overrides=refold_af3_overrides,
                    refold_b2_overrides=refold_b2_overrides,
                    refold_esm_overrides=refold_esm_overrides,
                )
            if res.ok:
                st.code(res.stdout or "(no output)", language="bash")
            else:
                st.error(f"launch.py exited {res.returncode}")
                st.code((res.stderr or res.stdout or "").strip() or "(no output)")

if launch_clicked and _ready_to_submit():
    results = []
    with st.spinner("Submitting to Euler..."):
        for job_label, job_config, overrides, job_num_runs in build_jobs(selected):
            res = euler.sbatch_via_launcher(
                remote_cwd, job_config, overrides, label, auto_refold, dry_run=False,
                num_runs=int(job_num_runs),
                refold_af3_overrides=refold_af3_overrides,
                refold_b2_overrides=refold_b2_overrides,
                refold_esm_overrides=refold_esm_overrides,
            )
            results.append((job_label, {
                "ok": res.ok,
                "jobid": res.jobid or "?",
                "rc": res.returncode,
                "refolds": res.refold_jobids,
                "err": (res.stderr or res.stdout or "").strip(),
                "detached": res.detached,
                "log_path": res.log_path,
            }))
    st.session_state["launch_results"] = results
    st.cache_data.clear()   # so the SLURM queue re-fetches with the new jobs
    st.rerun()
