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
| `hw_probe_rate.m` | beam-training rate sweeps (ideal serving beam — see gap below) |
| `ddbs_beam_hybrid.m` | **proposed** hybrid TD-PS split (PS absorbs the fc term) |
| `hw_probe_hybrid.m` | naive vs hybrid focusing gain across TD bits |
| `ttd_beam_impaired.m` | serving beam built with the same impaired hardware |
| `training_ddbs_e2e.m` | DDBS training with **both** stages impaired |
| `hw_probe_rate_e2e.m` | **money figure**: naive vs hybrid rate, both stages impaired |

## Report back
The TD-bit and σ_τ rate curves from `hw_probe_rate` (how far rate falls, and where)
plus the gain table — send those and I'll help formalize the robust-design lever.

---

## UPDATE — user's measured rate results + the verified fix

### Confirmed on the paper's own metric (SNR=15 dB, K=3, single-path)
ideal DDBS = 4.764, Perfect CSI = 5.028 bit/s/Hz.

| impairment | rate | % of ideal |
|---|---|---|
| B_td = 12 / 10 / 8 | 0.284 / 0.247 / 0.268 | **~6%** (random-pointing floor) |
| sigma_tau = 1 / 2 / 5 ps | 4.746 / 4.739 / 4.684 | 100 / 99 / 98% |
| B_ps = 3 / 2 | 4.772 / 4.738 | 100 / 99% |

The asymmetry holds end-to-end: **TD resolution is the single fragile axis.**

### Modeling gap found and fixed
`hw_probe_rate.m` uses the baseline `training_near_rainbow_2d`, whose serving beam
`TTD_beam` is generated **ideally**. So it measures only the estimation stage —
which is an argmax, hence robust to impairments that degrade the beam but preserve
the peak subcarrier (why jitter reads 98% on rate while focusing gain falls to
0.63). `training_ddbs_e2e.m` + `ttd_beam_impaired.m` now impair **both** stages
with one shared hardware delay LSB; use `hw_probe_rate_e2e.m` for the real figure.

Related asymmetry worth stating in the paper: a serving beam focuses at
sin(theta) in [-1,1] so its delay range is ~2 ns, while a DDBS **training** beam
carries the large intercept theta_t (~ -33) and spans ~140 ns. At equal hardware
LSB the training waveform is the fragile one — the problem is specific to DDBS
training, not near-field focusing in general.

### The fix, verified (`ddbs_beam_hybrid.m`, `hw_probe_hybrid.m`)
Split the per-antenna phase at the carrier:

    2*pi*f_m*tau_n  =  2*pi*fc*tau_n  +  2*pi*(f_m - fc)*tau_n

The first term is frequency-independent -> realized **exactly by the phase
shifter** (mod 2pi). The TTD then only carries the residual, so a delay error
d_tau produces at most pi*B*d_tau of phase error instead of ~2*pi*fc*d_tau:
**TD precision now scales with bandwidth B, not carrier fc** (2*fc/B ~ 12x).

Mean focus gain (ideal 0.975):

| B_td | 16 | 14 | 13 | 12 | 11 | 10 |
|---|---|---|---|---|---|---|
| naive | 0.968 | 0.877 | 0.581 | 0.096 | 0.070 | 0.083 |
| **hybrid** | 0.975 | 0.975 | 0.974 | **0.971** | **0.960** | 0.917 |

**~4.5 bits saved (~20x resolution relaxation)**, and a 3-bit PS still yields
~0.98 of the hybrid gain — the burden is **not** shifted onto the phase shifter.

### Next
Run `hw_probe_hybrid` (reproduce the table) and `hw_probe_rate_e2e` (rate version,
both stages impaired). Then formalize: analytic phase-error bound and the required
bits as a function of (Nt, B, fc, theta_t), plus the TD-count/precision trade
against the AoSA baseline.

