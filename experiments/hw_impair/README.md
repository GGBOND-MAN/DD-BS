# (B) Hardware-impairment robustness of DDBS — feasibility probe (MATLAB)

**Improve-on-existing target.** OJ-COMS 2026 (Qaid, Nasir, Al-Ahmadi, Y. Liu) reduced
the **TD count** with an AoSA/grouped-TTD architecture but assumed **ideal TD/PS
values**. The orthogonal, untouched axis is **TD/PS resolution & error** — what
finite-bit hardware and fabrication tolerances do to the DDBS multi-strip pattern.
This probe measures that, and points straight at the contribution.

## Pre-validated finding (full Nt=256, M=1024 — run to reproduce)

Mean realized focusing gain (ideal ≈ 0.97):

| impairment | result | read |
|---|---|---|
| **TD bits** (quantise the 140 ns delay range) | 8-bit → **0.07** (destroyed); needs ~13-bit | **catastrophic** |
| **TD jitter** σ_τ (fabrication, hits fixed lines) | <1 dB needs σ_τ ≲ 2 ps; 5 ps → 0.63; 10 ps → 0.18 | **tight** |
| **PS bits** | 2-bit → 0.88; ≥3-bit → ≥0.95 | **robust** |
| **phase noise** | 20° → 0.92; 30° → 0.85 | moderate |

**The fragility is asymmetric and concentrated in TD precision.** Root cause: the
DDBS delay range (~140 ns, dominated by the large linear intercept n·d·θ_t) is
thousands of carrier periods (fH period ≈ 31 ps), so low-resolution/erroneous TD
scrambles the phase. This is exactly why the baseline resorts to *fixed* delay
lines — and it is the gap a robust design should close.

**The contribution (concrete):** make DDBS tolerate low-resolution TD —
- **delay-range reduction**: fold the large linear intercept θ_t into the PS /
  exploit the angle periodicity so the residual delay range fits a modest-bit TTD;
- **coarse-PS + fine-TTD** split, or **sub-array TTD**;
- then quantify the TD-bit / ps-precision saving vs. rate, benchmarked against
  ideal DDBS and the AoSA baseline (the paper's own rate metric).

## Setup
Needs the extracted baseline (`cal_loc`, `delay_polar_2d`, `near_field_channel`,
`training_near_rainbow_2d`) at `<repo>/baseline_code/code_nf_distance_dependent_rainbow/`.
Edit `BASELINE_DIR` at the top of each script otherwise.

## Run
```matlab
cd experiments/hw_impair
hw_probe_gain     % focusing gain vs TD bits / TD jitter / PS bits / phase noise
hw_probe_rate     % beam-training average rate vs the same (paper's metric)
```

## Files
| file | role |
|---|---|
| `ddbs_beam_impaired.m` | DDBS beams with finite-bit TD/PS, TD jitter, phase noise, amp error (ideal ⇒ = baseline `delay_polar_2d`) |
| `hw_probe_gain.m` | focusing-gain sensitivity sweeps |
| `hw_probe_rate.m` | beam-training rate sweeps (reuses baseline training) |

## Report back
The TD-bit and σ_τ rate curves from `hw_probe_rate` (how far rate falls, and where)
plus the gain table — send those and I'll help formalize the robust-design lever.
