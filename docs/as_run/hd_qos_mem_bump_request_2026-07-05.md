# Request: raise es_wicky gpuhe QOS per-user memory cap (unblock idle GPUs)

**To:** Euler cluster support / es_wicky QOS coordinator
**From:** lvalsecchi (Leo Valsecchi), Wicky lab
**Date:** 2026-07-05

## The ask (two one-line changes)

Please raise the **per-user memory** ceiling (`MaxTRESPerUser`) on two es_wicky GPU QOS so my
jobs can actually use the GPU slots the same QOS already grants me. **GPU and CPU limits
unchanged** — memory only:

```
sacctmgr modify qos es_wicky/gpuhe/24  set MaxTRESPerUser=cpu=323,gres/gpu=50,mem=1400G
sacctmgr modify qos es_wicky/gpuhe/120 set MaxTRESPerUser=cpu=184,gres/gpu=29,mem=850G
```

(current values: `.../24` = `cpu=323,gres/gpu=50,mem=929G`; `.../120` = `cpu=184,gres/gpu=29,mem=531G`)

## Why

I'm running a GPU design campaign (HalluDesign generation) on the wicky gpuhe pool. Each job
uses **1 GPU + ~22–27 GB RAM** (measured peak MaxRSS: median 22.2G, p90 26.1G, max 27.2G;
requested 28G). The **memory** ceiling is binding **before** the GPU ceiling:

| QOS | GPU cap | mem cap | jobs the mem cap allows (@28G) | GPUs left idle |
|---|---|---|---|---|
| es_wicky/gpuhe/24 | 50 | 929G | **33** | **17** |
| es_wicky/gpuhe/120 | 29 | 531G | **18** | **11** |

So ~28 GPUs I'm already allotted sit unusable because the per-user mem cap tops out first
(all pending jobs are blocked on `QOSMaxMemoryPerUser`, not fairshare — my fairshare is 1.0,
RawUsage 0). Setting mem = GPU_cap × 28G (1400G / 812G→850G) lets the GPU cap be the real limit,
roughly **51 → 79 concurrent** for me, with **no change to the GPU/CPU allocation**.

This is for the current campaign; happy to have it reverted afterward if preferred.

Thanks!

---

### Notes for Leo (not part of the message)
- This is per-*user* (`MaxTRESPerUser`), so it only affects your own concurrency, not the pool total.
- If they'd rather not raise it, the fallback is time (HD finishes at ~51 concurrent + gpupr
  overflow in ~a day) or a gpuhe.bulk pilot (preemptible — needs testing first).
- Verify after they apply: `sacctmgr -n -P show qos es_wicky/gpuhe/24 format=Name,MaxTRESPU`
