# Job Distributor — concept draft

A periodic loop that rebalances **pending** jobs across the gpuhe queues so that
generators always win the high-priority tiers and refolds soak up whatever is
left — maximising throughput without starving the gen→refold dependency chain.

Status: **concept for review** (not built). Open questions at the end.

---

## 0. Core insight (why this works the way it does)
`gpuhe.4h / .24h / .120h / .bulk` are **not separate hardware** — all four are
the *same 145 GPU nodes*, differing only in:
- **PriorityTier**: 4h/24h/120h = **1000** (high), bulk = **250** (low, preemptible)
- **MaxTime**: 4h=4h, 24h=2d, 120h=15d, bulk=15d
- **per-tier QOS GPU ceiling** (es_wicky): /4=96, /24=50, /120=29, /bulk=240

So "distributing across queues" is really **priority-lane routing**: which lane a
job sits in decides how hard it competes for the shared GPUs. The distributor's
job is to keep generators in the fast lanes and park refolds in slower lanes
while gens are pending, then open all lanes to refolds once gens drain.

## 1. Rules (from the spec)
1. **No generators on bulk.** Gens only ever sit in 4h/24h/120h.
2. **Generators are prioritised** in 4h/24h/120h.
3. **While any generator is pending**, refolds go to **bulk only**.
   **Once no generator is pending**, refolds may use **all** queues (4/24/120h + bulk).
4. **Walltime tracks the queue**: 4h→`04:00:00`, 24h→`24:00:00`, 120h→`120:00:00`
   (bulk→`120:00:00` / max).

## 2. Hard constraints (discovered this session — the design must respect these)
- **Submit-time remap**: `--partition=gpuhe.bulk` (full name) is rewritten to
  `gpuhe.24h`; `--partition=bulk` (short) lands on bulk. → the distributor must
  **act post-submit via `scontrol`**, not rely on submit strings.
- **`scontrol update partition=` STICKS** (post-submit bypasses the remap plugin).
  Verified: moved 1 211 esm jobs onto bulk this way.
- **`scontrol` can only *decrease* `TimeLimit`** for non-admins (increase fails).
  ⇒ **submit every distributable job with a generous walltime (≥120h)** so the
  loop can place it on *any* tier by *reducing* TimeLimit to that tier's max.
  (A job submitted at 12h can never be moved up to the 120h walltime.)
- **GPU eligibility**: ESM + AF3 need **Ampere+ (sm_80+, eu-g4/g6 nodelist)**;
  B2 + all generators run on any 24G card. The loop must keep each job's nodelist/
  gres and only place it where it can actually run.
- **Coexistence with the nice-looper** (`refold_priority_babysitter`): that loop
  nices *bulk* refolds by progress. The distributor decides *which queue*; the
  babysitter decides *order within bulk*. Keep them separate (distributor must not
  clobber nice values; babysitter must not move partitions).

## 3. Architecture
A single long-running loop (mirror `refold_prio_loop.slurm` / the babysitter),
dry-run by default, `--go` to act, ~60–120 s cadence:

```
loop:
  jobs   = scontrol show job -o   (mine, PENDING)         # one cheap snapshot
  classify each: kind(gen|refold), refolder(af3|b2|esm),
                 cur_partition, gpu_need(ampere|any), walltime
  gens_pending = any(kind==gen and PENDING)
  empt[T]      = emptiness(T) for T in {4h,24h,120h,bulk}  # §4
  for each pending job (priority order §5):
     target = choose_queue(job, gens_pending, empt)        # §5 rules
     if target != cur_partition or timelimit != maxtime(target):
        scontrol update jobid=J partition=target TimeLimit=maxtime(target)
        (only if target is eligible for the job's gpu_need)
```

## 4. "Emptiness" metric (pluggable)
Since tiers share nodes, emptiness = *how readily a job placed there starts*.
Candidate signals (start simple, tune):
- `running_in[T] / (running_in[T] + pending_in[T])` for my jobs — high ⇒ empty.
- headroom vs the per-tier QOS GPU ceiling (/4=96,/24=50,/120=29,/bulk=240).
- recent start-rate per tier (jobs that left PENDING in the last cycle).
MVP: rank tiers by `(QOS_ceiling[T] − my_running[T])`, i.e. how much of my
per-tier quota is unused — fill the emptiest first.

## 5. Placement decision (`choose_queue`)
Priority order to process jobs: **generators first**, then refolds.
```
if job.kind == gen:
    # rule 1+2: never bulk; pick emptiest of {4h,24h,120h} the job FITS
    #   (fits = queue MaxTime >= job's required walltime; e.g. bd1~7h -> not 4h)
    return emptiest_fitting(job, [4h,24h,120h])
else:  # refold
    if gens_pending:
        return bulk                      # rule 3a: refolds to bulk while gens wait
    else:
        return emptiest_fitting(job, [4h,24h,120h,bulk])   # rule 3b: open all lanes
```
`emptiest_fitting` also enforces **gpu_need** (ESM/AF3 → only Ampere-backed
placements; respect the existing nodelist) and sets `TimeLimit = maxtime(T)`.

## 6. Walltime handling
- maxtime(4h)=04:00:00, maxtime(24h)=24:00:00, maxtime(120h)=120:00:00,
  maxtime(bulk)=120:00:00 (or 15d).
- Because TimeLimit can only shrink: **submit policy change** — have launch.py
  submit gens + refolds with a long ceiling (e.g. 120:00:00) so the distributor
  can place them anywhere. The distributor then *reduces* to the tier max.
- A job whose real runtime exceeds a tier's MaxTime is simply never placed there
  (e.g. a 7 h generator can't go to 4h).

## 7. Safety / anti-thrash
- **Hysteresis**: only move a job if the target tier is meaningfully emptier
  (e.g. ≥ threshold) — avoid ping-ponging each cycle.
- **Cap moves/cycle** (e.g. ≤200 `scontrol` calls) to spare slurmctld.
- **Never touch RUNNING jobs** (only PENDING) — moving a running job is disruptive.
- **Idempotent**: a job already on its target tier with the right TimeLimit is skipped.
- **Don't fight the babysitter**: distributor sets partition+TimeLimit only; the
  babysitter owns `nice`. (If both touch bulk refolds, that's fine — different fields.)
- **Dry-run default + log every move**; `--go` to act.

## 8. Edge cases / open questions
1. **Emptiness metric**: which signal do you want (QOS-headroom / running-ratio /
   start-rate)? MVP = QOS-headroom.
2. **bulk walltime**: 120h or the full 15d?
3. **Submit-walltime change**: OK to switch launch.py to submit gens+refolds at a
   120h ceiling (needed so the loop can place on any tier)? Without it, the loop
   can only move jobs *down*-tier, not up.
4. **Gens that exceed 24h** (bd1 ~7h is fine; anything >2d?) — confirm no gen needs
   >120h.
5. **Scope**: refolds = af3+b2+esm all distributed? (af3 needs Ampere; b2 any.)
6. **Cadence + move cap**: 60 s / 200 moves OK?
7. **Relationship to the existing 2 loops** (nice-looper + refold prio): fold the
   distributor INTO the babysitter loop, or run as a separate `job_distributor.sh`
   loop? (Separate is cleaner; both read the same scontrol snapshot.)

## 9. Phased build (once concept is approved)
- **P1**: read-only reporter — classify jobs, compute emptiness, print the moves it
  *would* make (validate the logic against live state).
- **P2**: generator protection — enforce rules 1–3 (gens off bulk, refolds→bulk while
  gens pending) via `scontrol partition=`, dry-run→`--go`.
- **P3**: emptiness-based spreading + walltime-per-tier (needs the submit-walltime change).
- **P4**: wrap in a babysitter-style loop + add to the existing prio-loop slurm job.
