"""
Per-stack metadata: which override fields exist, defaults, helptext.

The "path" field is the dotted YAML path that launch.py --override expects.
Codifies the per-stack syntax we hard-won; see memory entry on launch.py overrides.
"""
from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class ParamSpec:
    name: str           # display label
    path: str           # dotted override path
    default: str        # default value as string (we pass strings to launch.py anyway)
    help: str = ""
    choices: list[str] | None = None   # if set, render a dropdown instead of free text
    choice_labels: dict[str, str] | None = None  # value -> friendly dropdown label
    is_bool: bool = False              # render a checkbox; emits "true"/"false"
    as_list: bool = False              # comma-separated input -> JSON list override
    disabled: bool = False             # render the widget greyed-out (non-editable)
    info_only: bool = False            # display-only: render the widget but emit NO override
    group: str = ""                    # "" = always shown; "slurm"/"limpnn" go in expanders
    paper: str = ""                    # repo/paper-default note shown in the info box;
                                       # the "Paper / repo defaults" preset (app.py) sets these


@dataclass
class StackSpec:
    key: str            # campaigns/<target>/<key>/ directory name
    label: str          # display label
    config: str         # path to default.yaml (relative to rna_binder/)
    params: list[ParamSpec] = field(default_factory=list)


# GENERATION stacks PIN a specific GPU type to exclude the cards that break our
# models: a bare `--gpus=1` / `--gres=gpumem:24g` will happily land on a CC 7.5
# Turing card (titan_rtx / quadro_rtx_6000, which live in the gpuhe pool too) or a
# Blackwell rtx_pro_6000 — both of which silently break AF3/Triton/Protenix (and
# Blackwell breaks the bundled PyTorch for every stack). These are exactly the
# CC < 8.0 (and CC 12.0) cards. Verified 2026-06-09. Format "<type>:<count>".
#
# REFOLDS instead use a node-range ALLOW-LIST (`--nodelist=<Ampere/Ada nodes>`),
# which excludes the same bad cards but keeps the whole 3090+4090 pool eligible
# (shorter queues) — see REFOLD_GPU_STRATEGIES + expand_refold_gpu() in app.py.
# (The same node-range approach could replace the type-pin here for generation.)
GPU_CHOICES: list[str] = [
    "rtx_4090:1",   # 24 GB, Ada    (CC 8.9)
    "rtx_3090:1",   # 24 GB, Ampere (CC 8.6)
    "ampere_24g",   # 24 GB Ampere/Ada pool: 3090+4090 ONLY (no Turing). pair_gpumem
                    # expands to gpus=1 + gpumem:20g + AMPERE_PLUS_NODES nodelist.
                    # REQUIRED for HalluDesign (Protenix Evoformer needs sm_80+; Turing fails).
    "turing_24g",   # 24 GB "all usable" pool: 3090+4090+titan+quadro. app.pair_gpumem
                    # expands to gpus=1 + gpumem floor + the ALL_24G_NODES nodelist, so
                    # SLURM grabs whichever 24 GB card is free first (shortest queue).
                    # Generators run on Turing sm_75 with no flag (GENERATOR_GPU_COMPAT.md).
                    # DEFAULT for rfd3/boltzgen/boltzdesign. NOT for HD (Turing breaks Protenix).
    "a100:1",       # 40 GB, Ampere (CC 8.0)
    "a100_80gb:1",  # 80 GB, Ampere (CC 8.0)
]

# VRAM-forward dropdown labels (value stays the pinned type the override needs).
GPU_LABELS: dict[str, str] = {
    "rtx_4090:1":  "24 GB · rtx_4090",
    "rtx_3090:1":  "24 GB · rtx_3090",
    "ampere_24g":  "24 GB · Ampere only (3090/4090)",
    "turing_24g":  "24 GB · all usable (3090/4090/titan/quadro)",
    "a100:1":      "40 GB · a100",
    "a100_80gb:1": "80 GB · a100 (80 GB)",
}

# Euler's `gpumem` gres is a MINIMUM-total-VRAM filter, so it must be <= the
# chosen card's memory or the allocation is impossible ("Requested node
# configuration is not available"). The GPU dropdown only sets slurm.gpus, so
# we pair it with a matching slurm.gres derived from this map.
GPU_GPUMEM: dict[str, str] = {
    "rtx_4090:1":  "gpumem:20g",   # 24 GB card
    "rtx_3090:1":  "gpumem:20g",   # 24 GB card
    "a100:1":      "gpumem:40g",   # 40 GB card
    "a100_80gb:1": "gpumem:80g",   # 80 GB card
}

# HalluDesign runs Protenix internally (DeepSpeed Evoformer needs sm_80+), so it
# MUST stay on Ampere/Ada — NEVER Turing. Restrict HD's GPU dropdown (no turing_24g);
# 24 GB 3090/4090 are fine (HD was never VRAM-bound, only architecture-bound), a100 too.
HD_GPU_CHOICES: list[str] = ["ampere_24g", "rtx_4090:1", "rtx_3090:1", "a100:1", "a100_80gb:1"]

# Walltime fully determines the gpuhe partition tier (Euler's lua picks
# gpuhe.4h / .24h / .120h from --time), so the time field is a 3-tier dropdown.
TIME_CHOICES: list[str] = ["04:00:00", "24:00:00", "120:00:00"]
TIME_LABELS: dict[str, str] = {
    "04:00:00":  "4 h  · gpuhe.4h (fastest queue)",
    "24:00:00":  "24 h · gpuhe.24h",
    "120:00:00": "120 h · gpuhe.120h",
}

# Euler's job_submit.lua auto-selects the gpuhe.{4h,24h,120h} tier from the
# requested --time and OVERRIDES any explicit --partition, so a partition
# dropdown is pointless — walltime alone decides the tier. The one real
# partition lever is the low-priority, preemptible "bulk" queue, entered with
# `-p bulk`: a lua alias that maps GPU jobs to gpuhe.bulk at ANY walltime.
# (`-p gpuhe.bulk` does NOT work — the plugin rewrites it back to a normal tier,
# which is why the old "gpuhe.bulk" override silently ran on gpuhe.4h/24h.)
# The "Run on bulk" checkboxes in app.py set slurm.partition=bulk; every sbatch
# template already carries --requeue, which bulk needs to survive preemption.
BULK_PARTITION = "bulk"


STACKS: list[StackSpec] = [
    StackSpec(
        # Points at the new multif config (rfd3__multif) so launches go through
        # the external-LiMPNN path with mpnn_cli knobs, not the bundled-LiMPNN
        # path. Campaign dir layout is unchanged (rfd3__ligandmpnn/<label>_<uuid>/).
        key="rfd3__ligandmpnn",
        label="RFD3 + LigandMPNN",
        config="code/configs/rfd3__multif/default.yaml",
        params=[
            ParamSpec("n_batches", "generator_params.n_batches", "1",
                      help="RFD3 batches. Backbones = n_batches × diffusion_batch_size. "
                           "For exactly N backbones set n_batches = ceil(N / diffusion_batch_size) "
                           "(e.g. with batch 8: 200 → 25 batches). "
                           "Total seqs = backbones × (LiMPNN batch_size × number_of_batches).",
                      paper="Repo default 1; pure throughput knob, no paper value."),
            ParamSpec("diffusion_batch_size", "generator_params.diffusion_batch_size", "8",
                      help="RFD3 designs denoised in parallel per batch (RFD3 default 8 — a "
                           "VRAM/throughput knob, not a paper value). Backbones = n_batches × this; "
                           "lower for tight VRAM, raise for throughput.",
                      paper="RFD3 config default: 8 (RosettaCommons foundry "
                            "rfdiffusion3.yaml; 200 diffusion timesteps)."),
            ParamSpec("LiMPNN batch_size", "generator_params.mpnn_cli.batch_size", "4",
                      help="LigandMPNN sequences per batch. Total seqs/backbone = "
                           "batch_size × number_of_batches.",
                      paper="RFdiffusion→MPNN convention ≈8 seqs/backbone "
                            "(Bennett 2023 used 2; LigandMPNN argparse default 1).",
                      group="limpnn"),
            ParamSpec("LiMPNN number_of_batches", "generator_params.mpnn_cli.number_of_batches", "1",
                      help="LigandMPNN batches per backbone.",
                      paper="LigandMPNN argparse default: 1.",
                      group="limpnn"),
            ParamSpec("LiMPNN temperature", "generator_params.mpnn_cli.temperature", "0.1",
                      help="LigandMPNN sampling temperature.",
                      paper="Repo default & field standard: 0.1 "
                            "(ProteinMPNN suggests 0.1–0.3 for diversity).",
                      group="limpnn"),
            ParamSpec("gpu", "slurm.gpus", "turing_24g", choices=GPU_CHOICES, choice_labels=GPU_LABELS, group="slurm"),
            ParamSpec("slurm time", "slurm.time", "24:00:00", choices=TIME_CHOICES, choice_labels=TIME_LABELS, group="slurm"),
        ],
    ),
    StackSpec(
        key="halludesign__ligandmpnn",
        label="HalluDesign + LigandMPNN",
        config="code/configs/halludesign/default.yaml",
        params=[
            ParamSpec("num_seqs", "generator_params.params.num_seqs", "10",
                      help="Sequences sampled per recycle cycle. Total designs "
                           "≈ num_seqs × num_recycles.",
                      paper="HalluDesign argparse default: 8 (num_seqs)."),
            ParamSpec("num_recycles", "generator_params.params.num_recycles", "20",
                      help="Number of optimization cycles. Total designs "
                           "≈ num_seqs × num_recycles.",
                      paper="HalluDesign argparse default: 10; design examples use 10–30 "
                            "(Paper-defaults preset uses 15)."),
            ParamSpec("ref_time_steps", "generator_params.params.ref_time_steps", "30",
                      help="Diffusion truncation steps. Lower = local refinement, "
                           "higher = broader exploration. For structured RNA targets, "
                           "lower values produce better RMSD.",
                      paper="Passed to Protenix as diffusion_steps: 50 = local optimization "
                            "(argparse default), 150 = global optimization & design (every "
                            "README design example uses 150), 200 = full/untruncated diffusion "
                            "('pure prediction'). Paper-defaults preset uses 60 (live config)."),
            ParamSpec("num_packs", "generator_params.params.mpnn_num_packs", "1",
                      help="LigandMPNN sidechain packs per design (number_of_packs_per_design). "
                           "Total designs ≈ num_seqs × num_recycles × num_packs.",
                      paper="HalluDesign argparse default: 1 (mpnn_num_packs)."),
            # HD runs Protenix internally (DeepSpeed Evoformer needs Ampere+/sm_80+);
            # on Turing it errors "Unsupported GPU and data type combination" and silently
            # emits RNA-only designs. Restrict to Ampere/Ada (HD_GPU_CHOICES, no Turing);
            # default to the 24 GB Ampere pool (3090/4090). pair_gpumem sets gres + nodelist.
            ParamSpec("gpu", "slurm.gpus", "ampere_24g", choices=HD_GPU_CHOICES, choice_labels=GPU_LABELS, group="slurm"),
            ParamSpec("slurm time", "slurm.time", "04:00:00", choices=TIME_CHOICES, choice_labels=TIME_LABELS, group="slurm"),
        ],
    ),
    StackSpec(
        key="boltzgen__boltzif",
        label="BoltzGen + BoltzIF/LiMPNN (multif)",
        config="code/configs/boltzgen__multif/default.yaml",
        params=[
            # multif config reads bg_cli, NOT cli.flags (which is the non-multif
            # boltzgen path). Verified against build_boltzgen_multif_sbatch.
            # BoltzGen multif is pair-only: launch.py requires inverse_folders
            # == {boltzif, ligandmpnn}, so this single stack always runs BOTH
            # inverse folders on one shared backbone — hence both the bif_cli
            # and mpnn_cli knobs below. Outputs land in boltzgen__boltzif/ AND
            # boltzgen__ligandmpnn/ (plus the shared _backbones/ dir).
            ParamSpec("num_designs (backbones per run)", "generator_params.bg_cli.num_designs", "1",
                      help="Number of protein backbones generated per run. "
                           "Total backbones = num_designs × num_runs. Production: 200.",
                      paper="BoltzGen repo default: 10000; README recommends 10k–60k for "
                            "production (test first at 50). Scale knob — left out of the preset."),
            ParamSpec("BIF seqs/backbone", "generator_params.bif_cli.inverse_fold_num_sequences", "4",
                      help="BoltzIF sequences per backbone.",
                      paper="BoltzGen repo default: 1 seq/backbone "
                            "(inverse_fold_num_sequences). BoltzIF exposes no temperature.",
                      group="bif"),
            ParamSpec("LiMPNN batch_size", "generator_params.mpnn_cli.batch_size", "4",
                      help="LigandMPNN sequences per batch. Total LiMPNN seqs/backbone = batch_size × number_of_batches.",
                      paper="RFdiffusion→MPNN convention ≈8 seqs/backbone "
                            "(Bennett 2023 used 2; LigandMPNN argparse default 1).",
                      group="limpnn"),
            ParamSpec("LiMPNN number_of_batches", "generator_params.mpnn_cli.number_of_batches", "1",
                      help="LigandMPNN batches per backbone.",
                      paper="LigandMPNN argparse default: 1.",
                      group="limpnn"),
            ParamSpec("LiMPNN temperature", "generator_params.mpnn_cli.temperature", "0.1",
                      help="LigandMPNN sampling temperature.",
                      paper="Repo default & field standard: 0.1 "
                            "(ProteinMPNN suggests 0.1–0.3 for diversity).",
                      group="limpnn"),
            ParamSpec("gpu", "slurm.gpus", "turing_24g", choices=GPU_CHOICES, choice_labels=GPU_LABELS, group="slurm"),
            ParamSpec("slurm time", "slurm.time", "04:00:00", choices=TIME_CHOICES, choice_labels=TIME_LABELS, group="slurm"),
        ],
    ),
    StackSpec(
        # Points at the new multif config — BD1 runs hallucination only
        # (--run_ligandmpnn false --run_alphafold false set in YAML), then
        # external LiMPNN with mpnn_cli knobs. ~1.5x faster than the bundled
        # stack and customizable. Campaign dir layout is unchanged.
        key="boltzdesign__ligandmpnn",
        label="BoltzDesign1 + LigandMPNN",
        config="code/configs/boltzdesign__multif/default.yaml",
        params=[
            ParamSpec("design_samples", "generator_params.design_samples", "50",
                      help="Number of hallucination trajectories. Production ≥50; "
                           "smoke ≥5 (fewer often produce 0 survivors of BD1's "
                           "internal confidence filter).",
                      paper="BD1 argparse default: 1; README examples use 2–5. Confidence "
                            "filters i_ptm≥0.5 & complex_pLDDT≥0.7 are BD1 repo defaults. "
                            "Scale knob — left out of the preset."),
            ParamSpec("LiMPNN batch_size", "generator_params.mpnn_cli.batch_size", "4",
                      help="LigandMPNN sequences per batch. Total seqs/backbone = "
                           "batch_size × number_of_batches.",
                      paper="RFdiffusion→MPNN convention ≈8 seqs/backbone "
                            "(Bennett 2023 used 2; LigandMPNN argparse default 1).",
                      group="limpnn"),
            ParamSpec("LiMPNN number_of_batches", "generator_params.mpnn_cli.number_of_batches", "1",
                      help="LigandMPNN batches per backbone.",
                      paper="LigandMPNN argparse default: 1.",
                      group="limpnn"),
            ParamSpec("LiMPNN temperature", "generator_params.mpnn_cli.temperature", "0.1",
                      help="LigandMPNN sampling temperature.",
                      paper="Repo default & field standard: 0.1 "
                            "(ProteinMPNN suggests 0.1–0.3 for diversity).",
                      group="limpnn"),
            ParamSpec("gpu", "slurm.gpus", "turing_24g", choices=GPU_CHOICES, choice_labels=GPU_LABELS, group="slurm"),
            ParamSpec("slurm time", "slurm.time", "24:00:00", choices=TIME_CHOICES, choice_labels=TIME_LABELS, group="slurm"),
        ],
    ),
]


def by_key(key: str) -> StackSpec | None:
    for s in STACKS:
        if s.key == key:
            return s
    return None


# ── refold settings ─────────────────────────────────────────────────────────
# One global panel feeds these as --refold-af3-override / --refold-b2-override
# to launch.py (only applied when Auto-refold is on). Paths match the keys in
# code/configs/refold__af3/default.yaml and refold__b2/default.yaml.

# Config files the refold ParamSpec paths target (relative to the rna_binder
# project root / eulermount). Used to seed the refold widgets from the live YAML,
# the way StackSpec.config does for the per-stack widgets.
REFOLD_AF3_CONFIG = "code/configs/refold__af3/default.yaml"
REFOLD_B2_CONFIG = "code/configs/refold__b2/default.yaml"
REFOLD_ESM_CONFIG = "code/configs/refold__esm/default.yaml"

# ── refold GPU selection (supervisor-style: VRAM floor + node-range) ───────────
# Each strategy key expands (app.expand_refold_gpu) to slurm.gpus=1 +
# slurm.gres=gpumem:.. + optional slurm.nodelist (+ slurm.partition for a100),
# matching refold__{af3,b2} default.yaml. AF3's *default* Triton FlashAttention
# needs CC >= sm_80, BUT the portable xla attention path runs AF3 on Turing too:
# refold_af3.sh auto-picks xla when compute_cap < 80 (validated 2026-06-11 on
# titan_rtx + quadro_rtx_6000 — confidences match the 3090 within AF3's sampling
# noise, ranking_score Δ≤0.04, and it ran faster). The all-usable pool
# (turing_24g = Ampere/Ada + Turing) is now the AF3 DEFAULT; pick ampere_24g to
# restrict to 3090/4090. Boltz-2 runs on any GPU. The a100 strategies (40/80 GB, for
# large complexes) live in the gpupr partition — NOT gpuhe — so they set
# slurm.partition to the gpupr tier matching the walltime; gpupr holds only a100s,
# so a gpumem floor alone selects the SKU (no nodelist needed). See GPU_REFOLD_COMPAT.md.
AMPERE_PLUS_NODES = "eu-g4-[001-032],eu-g6-[001-080]"   # rtx_3090 + rtx_4090
RTX4090_NODES = "eu-g6-[001-080]"                        # rtx_4090 only
# Turing 24 GB cards (sm_75, AF3 via xla): quadro_rtx_6000 eu-g3-[062-066]/[083-084]
# + eu-lo-g3-[049-061]; titan_rtx eu-g3-[067-079] (062-079 is a contiguous quadro+titan
# block). Excludes the 11 GB 2080 Ti (eu-lo-g2) and Blackwell rtx_pro_6000 (eu-g7).
TURING_24G_NODES = "eu-g3-[062-079],eu-g3-[083-084],eu-lo-g3-[049-061]"
ALL_24G_NODES = AMPERE_PLUS_NODES + "," + TURING_24G_NODES   # Ampere/Ada + Turing

# key -> (label, gpumem gres, nodelist, partition)
#   partition "" = leave default (gpuhe, walltime-tiered); "gpupr" = a100 pool,
#   expanded to gpupr.<tier> from the selected walltime.
REFOLD_GPU_STRATEGIES: dict[str, tuple[str, str, str, str]] = {
    # Public gpu.* pool (Turing: 2080 Ti 11 GB + titan/quadro 24 GB) on a SEPARATE
    # fairshare lane from wicky gpuhe. The gpumem:8g floor is what KEEPS the job on
    # the public pool — a >=24 GB request is routed by Euler's plugin to wicky
    # gpuhe (verified 2026-06-19). partition "gpu" => gpu.<tier-by-walltime>.
    # Boltz-2 only (~2.7 GB needed); AF3/ESM need Ampere so they can't use this pool.
    "public_turing": ("Public gpu.* · Turing (Boltz-2, 8 GB floor)", "gpumem:8g", "", "gpu"),
    "ampere_24g":  ("24 GB · Ampere+ (3090/4090)", "gpumem:24g", AMPERE_PLUS_NODES, ""),
    "rtx4090_24g": ("24 GB · 4090 only",           "gpumem:24g", RTX4090_NODES, ""),
    "turing_24g":  ("24 GB · all usable (3090/4090/titan/quadro)", "gpumem:24g", ALL_24G_NODES, ""),
    "any_24g":     ("24 GB · any GPU",             "gpumem:24g", "", ""),
    "a100_40g":    ("40 GB · a100 (gpupr)",        "gpumem:40g", "", "gpupr"),
    "a100_80g":    ("80 GB · a100 80 GB (gpupr)",  "gpumem:80g", "", "gpupr"),
}
# AF3 default = ampere_24g (3090/4090 only). AF3 needs Ampere+ (Triton/Evoformer);
# the xla Turing path is NOT working in prod — AF3 fails exit-127 on titan/quadro/eu-g3
# (verified 2026-06-16). turing_24g stays selectable but will fail AF3. a100 also OK.
REFOLD_AF3_GPU_CHOICES = ["ampere_24g", "rtx4090_24g", "turing_24g", "a100_40g", "a100_80g"]
REFOLD_B2_GPU_CHOICES = ["public_turing", "any_24g", "ampere_24g", "rtx4090_24g", "a100_40g", "a100_80g"]
# ESMFold2's Triton pair-bias kernel needs sm_80+ (Ampere/Ada) — Turing fails, so
# no turing_24g/any_24g here. a100 (sm_80) works.
REFOLD_ESM_GPU_CHOICES = ["ampere_24g", "rtx4090_24g", "a100_40g", "a100_80g"]
REFOLD_GPU_LABELS = {k: v[0] for k, v in REFOLD_GPU_STRATEGIES.items()}

REFOLD_AF3_PARAMS: list[ParamSpec] = [
    ParamSpec("gpu", "slurm.gpus", "ampere_24g",
              choices=REFOLD_AF3_GPU_CHOICES, choice_labels=REFOLD_GPU_LABELS,
              help="Selects by VRAM + node range. Default = ampere_24g (3090/4090 only). "
                   "AF3 needs Ampere+ (Triton/Evoformer); on Turing (titan/quadro/eu-g3) it "
                   "fails exit-127 — the xla Turing path is NOT working in prod (verified "
                   "2026-06-16). a100 also OK. See GPU_REFOLD_COMPAT.md."),
    ParamSpec("slurm time", "slurm.time", "24:00:00", choices=TIME_CHOICES, choice_labels=TIME_LABELS),
    ParamSpec("model_seeds", "refolder_params.model_seeds", "1",
              help="Comma-separated AF3 seeds, e.g. 1,2,3", as_list=True,
              paper="AF3 uses the seeds you supply; paper recommends 5 seeds × 5 samples "
                    "= 25 predictions for confident calls. Pipeline uses 1 (one refold/design)."),
    ParamSpec("num_diffusion_samples", "refolder_params.num_diffusion_samples", "1",
              help="AF3 diffusion samples per design (1 = one refold per cif)",
              paper="AF3 repo default: 5 per seed. Pipeline deliberately uses 1 "
                    "(one refold per design) — not changed by the preset."),
    # AF3 always runs BOTH passes per design (launch.py sets msa_variants =
    # [False, True] and overwrites use_rna_msa per variant), so this is not a user
    # toggle. Show it as a fixed, disabled "both" indicator and emit no override.
    ParamSpec("RNA MSA", "refolder_params.use_rna_msa", "msaoff",
              choices=["msaoff"], choice_labels={"msaoff": "msaoff only (msaon deactivated)"},
              disabled=True, info_only=True,
              help="AF3 msaon (RNA-MSA) refolds are DEACTIVATED — only the msaoff "
                   "pass runs. Flip AF3_RUN_MSAON in launch.py to reactivate msaon."),
]

REFOLD_B2_PARAMS: list[ParamSpec] = [
    ParamSpec("gpu", "slurm.gpus", "public_turing",
              choices=REFOLD_B2_GPU_CHOICES, choice_labels=REFOLD_GPU_LABELS,
              help="Boltz-2 runs on any GPU (~2.7 GB). Default = public_turing: the "
                   "public gpu.* pool (idle 2080 Ti + titan/quadro) on a separate "
                   "fairshare lane, freeing wicky gpuhe for AF3/ESM. Pick any_24g/ampere "
                   "to put it back on gpuhe."),
    ParamSpec("slurm time", "slurm.time", "24:00:00", choices=TIME_CHOICES, choice_labels=TIME_LABELS),
    ParamSpec("diffusion_samples", "refolder_params.diffusion_samples", "1",
              help="Boltz-2 diffusion samples per design (1 = one refold per cif)",
              paper="Boltz-2 repo default: 1 (recycling 3, sampling 200). README suggests "
                    "recycling 10 + diffusion 25 for AF3-level accuracy. Pipeline uses 1."),
    # Boltz-2 has no RNA-MSA mode (launch.py forces its label to msaoff and never
    # applies use_rna_msa), so keep the control visible but disabled to make that
    # clear. Output-format and use_rna_template options removed — neither is usable
    # here; refolds always use the config default.
    ParamSpec("use_rna_msa", "refolder_params.use_rna_msa", "false",
              help="Boltz-2 has no RNA-MSA mode — not applicable.",
              is_bool=True, disabled=True),
]

# ── ESMFold2 (refold_esm) — local biohub/ESMFold2, single-sequence (msaoff) ─────
REFOLD_ESM_PARAMS: list[ParamSpec] = [
    ParamSpec("gpu", "slurm.gpus", "ampere_24g",
              choices=REFOLD_ESM_GPU_CHOICES, choice_labels=REFOLD_GPU_LABELS,
              help="ESMFold2's fused pair-bias attention is a Triton kernel needing "
                   "sm_80+ (Ampere/Ada); Turing fails. Default ampere_24g (3090/4090); "
                   "a100 also OK. ESMC-6B backbone uses ~21 GB."),
    ParamSpec("slurm time", "slurm.time", "24:00:00", choices=TIME_CHOICES, choice_labels=TIME_LABELS),
    ParamSpec("num_loops", "refolder_params.num_loops", "20",
              help="ESMFold2 folding loops (default 20)."),
    ParamSpec("num_sampling_steps", "refolder_params.num_sampling_steps", "100",
              help="ESMFold2 diffusion sampling steps (default 100). ~4 s/design at 100; "
                   "lower = faster, less precise."),
    ParamSpec("num_diffusion_samples", "refolder_params.num_diffusion_samples", "1",
              help="ESMFold2 samples per design (1 = one refold per cif)."),
    # ESM is single-sequence (msaoff) like B2; keep the control visible but disabled.
    ParamSpec("use_rna_msa", "refolder_params.use_rna_msa", "false",
              help="ESMFold2 refolds run single-sequence (msaoff). MSA path not wired.",
              is_bool=True, disabled=True),
]


def estimate_designs(stack_key: str, ov: dict[str, str]) -> tuple[int | None, str]:
    """Rough design-count estimate for the pre-launch scale readout (per run).

    ov maps ParamSpec.path -> current value. Returns (count, formula_string), or
    (None, "") for stacks we have no formula for. These are ESTIMATES — e.g. BIF
    caps at ~1 effective output and the HD count is approximate — so the GUI
    should label them as such and show the formula for sanity-checking.
    """
    def g(path: str, default: int) -> int:
        try:
            return int(float(ov.get(path, default)))
        except (TypeError, ValueError):
            return default

    if stack_key == "rfd3__ligandmpnn":
        nb = g("generator_params.n_batches", 1)
        dbs = g("generator_params.diffusion_batch_size", 8)
        bs = g("generator_params.mpnn_cli.batch_size", 4)
        nbt = g("generator_params.mpnn_cli.number_of_batches", 1)
        bb = nb * dbs
        return bb * bs * nbt, f"{nb} batches × {dbs} = {bb} backbones × {bs}×{nbt} seqs"
    if stack_key == "halludesign__ligandmpnn":
        ns = g("generator_params.params.num_seqs", 10)
        nr = g("generator_params.params.num_recycles", 20)
        npk = g("generator_params.params.mpnn_num_packs", 1)
        return ns * nr * npk, f"~{ns}×{nr}×{npk}"
    if stack_key == "boltzgen__boltzif":
        bb = g("generator_params.bg_cli.num_designs", 1)
        bif = g("generator_params.bif_cli.inverse_fold_num_sequences", 4)
        bs = g("generator_params.mpnn_cli.batch_size", 4)
        nbt = g("generator_params.mpnn_cli.number_of_batches", 1)
        return bb * (bif + bs * nbt), f"{bb}×({bif} BIF + {bs}×{nbt} LiMPNN)"
    if stack_key == "boltzdesign__ligandmpnn":
        ds = g("generator_params.design_samples", 200)
        bs = g("generator_params.mpnn_cli.batch_size", 4)
        nbt = g("generator_params.mpnn_cli.number_of_batches", 1)
        return ds * bs * nbt, f"{ds}×{bs}×{nbt}"
    return None, ""
