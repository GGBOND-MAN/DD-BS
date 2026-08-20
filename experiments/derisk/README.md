# Phase 0 — Multipath de-risk probe (MATLAB)

Purpose: **fail fast**. Before investing in the theory of the proposed direction
(*multipath-robust near-field beam training that exploits the distance-dependent
beam-split structure*), these scripts test the two empirical premises the whole
idea rests on. They reuse the baseline code unmodified and only swap the
single-path channel for an L-path one.

Baseline: T. Zheng, M. Cui, Z. Wu, L. Dai, *"Near-field wideband beam training
based on distance-dependent beam split,"* IEEE TWC 24(2), 2025.

---

## Setup

1. Extract the baseline code `code_nf_distance_dependent_rainbow.zip` (on the
   `main` branch of this repo) so that this folder sees it at:

   ```
   <repo>/baseline_code/code_nf_distance_dependent_rainbow/*.m
   <repo>/experiments/derisk/*.m        <- these scripts
   ```

   If you put the baseline elsewhere, just edit the `BASELINE_DIR` line at the
   top of `derisk_gateA.m` and `derisk_gateB.m`.

2. Requirements: MATLAB (R2020a, as the baseline), **Communications Toolbox**
   (`awgn`) and **Parallel Computing Toolbox** (`parfor`; if you don't have it,
   `parfor` still runs serially).

---

## Run

```matlab
cd experiments/derisk
derisk_gateA     % Gate A: does the baseline break under multipath?
derisk_gateB     % Gate B: does the DDBS structure make paths separable?
```

Each script saves a `*_results.mat` and pops the figures described below.
Default settings are tuned for a quick pass; scale up `N_iter` / dictionary
grids (`g1d`,`g2d`) once the qualitative answer is clear.

Rough runtime (Nt=256, M=1024, defaults): Gate A a few minutes (dominated by the
match-filter dictionary build + per-iteration MF correlation); Gate B similar.
Memory: Gate B's fine dictionary (`g1d=256, g2d=40`) is ~250 MB — lower `g1d` if
tight.

---

## What each gate decides

### Gate A — is the multipath weakness real? (`derisk_gateA.m`)
- **(A1)** Average rate vs number of paths `L` (fixed NLoS/LoS ratio, fixed SNR).
- **(A2)** Average rate vs NLoS/LoS power ratio (fixed `L`, fixed SNR).
- Curves: Perfect CSI (per-subcarrier analog MRT), DDBS on-grid, DDBS + match filter.

**GO** if DDBS on-grid and DDBS+MF fall clearly below Perfect CSI as `L` grows /
as the NLoS ratio rises. **NO-GO / re-scope** if the baseline stays near-optimal
even with strong NLoS (then the weakness isn't consequential in this regime —
pivot to a blockage / multi-beam-combining scenario, or to UPA).

Note: two distinct effects can open the gap and both are legitimate motivation —
(i) the single-location match filter *mis-locks* under template mismatch, and
(ii) even when it locks to the LoS path, a single beam cannot harvest the
multipath combining gain that Perfect CSI gets.

### Gate B — is the DDBS structure separable? (`derisk_gateB.m`)
- **(B1)** Match-filter correlation surface for one well-separated `L`-path
  realisation, with the true path locations overlaid. Look for `L` distinct
  maxima near the true paths.
- **(B2)** A deliberately-dumb greedy peak extractor on that surface; reports
  detection probability and localisation RMSE (angle + distance) vs SNR.

**GO** if ≥2 paths are reliably separated at moderate SNR with `K=3` pilots →
the "path = structured trajectory, resolvable by beam-split" premise holds →
proceed to design the real (structured / super-resolution) estimator. **Re-think**
if the surface fuses paths into one blob even at high SNR.

The greedy extractor is **not** the proposed algorithm — it only checks whether
the information needed to separate paths is present in the DDBS observation.

---

## Files
| file | role |
|---|---|
| `derisk_params.m` | system parameters (verbatim from `Rate_snr.m`) |
| `near_field_channel_multipath.m` | L-path near-field wideband channel (Eq. (2)) |
| `gen_multipath_scenario.m` | random LoS+NLoS geometry & gains |
| `build_match_dictionary.m` | single-location MF dictionary (as in `Rate_snr.m`) |
| `derisk_gateA.m` | Gate A driver + plots |
| `derisk_gateB.m` | Gate B driver + plots |

## Reporting back
Two numbers settle the go/no-go: **(A2)** the rate gap (Perfect − DDBS+MF) at
NLoS/LoS = −3 dB, and **(B2)** detection probability at SNR = 10 dB. Send those
two plus the (B1) surface and I'll help decide theory design vs re-scope.
