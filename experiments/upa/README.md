# UPA / 3D-DDBS — minimal feasibility probe (MATLAB)

Phase-0 for the pivot direction: extend distance-dependent beam split (DDBS) from
a ULA (1 angle + distance) to a **UPA / planar array** (azimuth + elevation +
distance). Collision-checked open lane; the novelty is the **distance/range
dimension + azimuth/elevation/range coupling** — the closest prior UPA work
(arXiv 2512.24727) is beam-squint **angle-only** sensing and never touches range.

These scripts de-risk the two UPA-specific questions before any theory.

## Setup
Needs the extracted baseline (for `cal_loc.m`, used by the coverage probe):
```
<repo>/baseline_code/code_nf_distance_dependent_rainbow/*.m
<repo>/experiments/upa/*.m
```
Edit `BASELINE_DIR` at the top of `upa_probe_angular_coverage.m` if placed elsewhere.
`upa_probe_coupling.m` is self-contained (no baseline needed).

## Run
```matlab
cd experiments/upa
upa_probe_coupling          % CORE: the near-field cross-coupling wall
upa_probe_angular_coverage  % SUPPORT: 2D angular coverage / pilot intuition
```

## The two questions

### CORE — can a CHEAP array focus in 3D? (`upa_probe_coupling.m`)
The UPA near-field path-length has a **non-separable cross term** −n_x n_y d²Ω_xΩ_y/r.
A cheap **separable TD architecture** (per-axis delays, N_x+N_y TDs) can match the
two angular linear terms and the two diagonal curvatures but **cannot realize the
cross term**; only a **fully-connected** TD array (N_x·N_y TDs, expensive) can. The
probe maps the focusing loss of dropping the cross term over the ±60° angular
region, at r parameterised by the **Rayleigh distance** (a full-Taylor curve, ≈0 dB
in the moderate near field, confirms the cross term is the sole culprit).

**Pre-validated finding (64×64, run it to reproduce + get the heatmaps):**

| r / Rayleigh | separable median loss | worst off-broadside corner | region within 3 dB |
|---|---|---|---|
| 0.05 (deep NF) | −1.4 dB | **−12 dB** | 67% |
| 0.12 | −0.2 dB | −2.8 dB | 100% |
| 0.30 (shallow NF) | −0.07 dB | −1.0 dB | 100% |

So it is **GREEN with a real, well-posed problem**, not a coin-flip:
- In the **moderate/shallow near field** the cheap **separable 3D-DDBS just works**
  (loss < 1 dB almost everywhere) → the direction is feasible.
- In the **deep near field (r ≲ 0.1·Rayleigh) at double-off-broadside angles**
  (both Ω_x and Ω_y large), the separable architecture **collapses (−12 to −16 dB)**
  — precisely where near-field focusing matters most. That corner is the paper's
  contribution: characterise the cross-coupling wall (it scales cleanly with
  r/Rayleigh and Ω_xΩ_y) and add a **low-cost correction** (a few joint/corner TDs,
  a semi-connected sub-array, or a coverage strategy that routes the deep-NF corner
  to a fully-connected sub-panel).

Read off the printed `worst` loss and `frac<3dB` per r, plus the heatmaps.

### SUPPORT — how many pilots to cover 2D angle? (`upa_probe_angular_coverage.m`)
On a UPA, one frequency axis makes a separable design trace a **1-D curve** through
the 2-D (Ω_x, Ω_y) plane, not an automatic fill. Covering it needs different
per-axis sweep rates (space-filling / Lissajous) + a few interleaved pilots. The
probe plots the locus and prints an approximate 2D coverage % for K pilots — the
UPA pilot-overhead intuition that is the whole point of a beam-split scheme.

## Report back
Two things settle the go/decide: **(core)** `frac<3dB` at r = 20 m and the loss
heatmaps; **(support)** the coverage % and locus plot. Send those and I'll help
pick the narrative and design the full scheme.

## Files
| file | role |
|---|---|
| `upa_delta_exact.m` | exact UPA near-field path-length difference |
| `upa_delta_taylor.m` | 2nd-order Taylor, with/without the cross term |
| `upa_probe_coupling.m` | CORE probe: separable-architecture focusing loss |
| `upa_probe_angular_coverage.m` | SUPPORT probe: 2D angular coverage locus |
