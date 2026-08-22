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

### CORRECTION — the first end-to-end probe was flawed

The user's first `hw_probe_rate_e2e` run showed **no** hybrid advantage
(B_td=12: naive 0.115, hybrid 0.166) despite the training-beam focusing gain
being 0.096 vs 0.971. Cause: the probe applied the fix to **training only** and
always served with a **pure-TD** `ttd_beam_impaired`, so the serving beam became
the bottleneck and dominated the rate.

Serving-beam focusing gain (its delay range is only ~2.1 ns):

| LSB [ps] | 16.5 | 32.9 | 131.7 |
|---|---|---|---|
| pure-TD | 0.650 | **0.044** | 0.043 |
| hybrid | 0.999 | **0.996** | 0.941 |

So the required **precision (LSB) is set by the carrier period, not by the delay
range** — even a 2.1 ns beam dies at a 32.9 ps LSB. The range only sets the
**bit count**, `bits = log2(range/LSB)`, which is why the 140 ns training
waveform needs ~6 more bits than the 2.1 ns serving beam at equal precision.

Fixed: `ttd_beam_impaired` takes a `use_hybrid` flag, `training_ddbs_e2e` passes
it through, and `hw_probe_rate_e2e` now runs **one architecture end-to-end per
arm** (naive trains+serves pure-TD; hybrid trains+serves with the TD-PS split).
Re-run `hw_probe_rate_e2e` for the corrected figure.

Also confirmed by the user: the PS does **not** become the new bottleneck —
at B_td=10, B_ps = Inf/6/5/4/3 gives 0.917/0.917/0.916/0.912/0.894.

### Next
Re-run `hw_probe_rate_e2e` (corrected). Then formalize: analytic phase-error bound
and required bits as a function of (Nt, B, fc, theta_t) — expect roughly
`bits >~ log2(Nt*d*|theta_t|*B/c)` for the hybrid vs the same with `fc` in place
of `B` for pure-TD — plus the TD-count x TD-precision cost surface against the
AoSA baseline.


---

## TASK 1 RESULT — redesigning DD-BS for a shared-TTD architecture (`redesign_for_shared.m`)

Rate (bit/s/Hz) with the recalibrated lookup, SNR=15 dB, Nt=256, M=1024:

| `s` | `K` | `theta_t` | range | 256 TTD | 128 TTD | 64 TTD | 32 TTD |
|---|---|---|---|---|---|---|---|
| 1.00 (baseline) | 3 | 31.73 | 134.8 ns | 4.740 | 4.689 | 3.665 | **1.683** |
| 0.50 | 6 | 15.86 | 67.4 ns | 4.788 | 4.771 | 4.686 | 4.127 |
| 0.25 | 12 | 7.93 | 33.7 ns | 4.744 | 4.733 | 4.721 | **4.767** |

**32 TTDs (8x fewer) + K=12 pilots reaches 100.6% of the full-per-antenna ideal**,
where the baseline parameter set gets 35%. Slowing the sweep costs no rate when
TTDs are plentiful (first column is flat), so this is a pure hardware-for-pilots
exchange. It also cuts the required delay range 4x (134.8 → 33.7 ns), attacking
the unrealizability problem of THEORY §13 with the same knob.

Frontier: `K_min ~ 1.5 P`, i.e. **`K * N_TTD ~ 1.5 Nt`**, holding at every entry.
Derivation and the falsifiable `K/P ∝ B/fc` prediction are in `THEORY.md` §14.

### Next check — `frontier_ttd_pilot.m`
Fills the frontier in finely (adds P=16) and **tests** the predicted scaling:
at P=8, `K_min` should be 6 / 12 / 24 for B = 2.5 / 5 / 10 GHz. Run it and send
the two tables; if the B-scaling holds, the law is a design rule rather than a
fitted constant.

---

## FRONTIER RESULT — `K_min = ceil(1.12 P)`, independent of bandwidth (`frontier_ttd_pilot.m`)

| P | N_TTD | K_min | `P*|theta_t|` |
|---|---|---|---|
| 2 | 128 | 3 (K>=3 floor) | 63.5 |
| 4 | 64 | 5 | 76.1 |
| 8 | 32 | 9 | 84.6 |
| 16 | **16** | 18 | 84.6 |

The invariant is `P*|theta_t| <= ~85`, giving `K_min = ceil(1.12 P)` — which
reproduces every measured point. **`K ~ 1.5P` in §14 was a coarse-ladder artifact.**

**The `B/fc` prediction was falsified**: at P=8, `K_min = 9` for B = 2.5 / 5 / 10 GHz
alike (predicted 6 / 12 / 24). The replacement mechanism is coverage truncation by
the tilted sub-array factor, in which both `fc` and `B` cancel; it predicts
`K_min = (0.84..1.04) P` from first principles vs the measured 1.12 P. Full
derivation in `THEORY.md` §15.

**Consequence:** `theta_t` and `theta_p` enter the design equation as a *ratio*, and
only `theta_t`'s DC re-centring role is large. So the pilot/TTD tradeoff is a
property of how DDBS places its sweep — potentially removable, not just tradeable.

### Next
| file | question |
|---|---|
| `decouple_theta.m` | Scale `theta_t`, `theta_p` independently at P=8, K=3. Does a small-`theta_t` / large-`theta_p` point reach ~100% of ideal **with full coverage**? If yes the tradeoff is broken. |
| `diag_shared_mechanism.m` | Measures usable bandwidth and coverage holes directly — confirms or refutes the truncation mechanism. |

---

## MECHANISM — coverage holes, and a deterministic design test

`diag_shared_mechanism.m` + `decouple_theta.m`:

- **Gain truncation refuted.** Measured usable bandwidth is 5-25x wider than the
  sub-array-factor prediction (2.64 vs 0.21 GHz at P=8,K=3).
- **Coverage holes confirmed, perfectly.** The largest interior gap of the
  recalibrated focus map, in array beamwidths `2/Nt`:
  `<= 1` → 99-100% of ideal (4/4); `>= 5` → 42-89% (5/5). No overlap, sharp cliff.
  `span` reads 1.04 even in total failure — it is useless; the **gap** is the statistic.
- **The `theta_t`/`theta_p` ratio equation is refuted.** At P=8, K=3 no independent
  `(theta_t, theta_p)` pair beats 46% of ideal, and shrinking `theta_t` alone makes it
  *worse* (12% at `s_t=1/16`). So `P*|theta_t| <= 85` is necessary but not sufficient —
  it holds only along the locus where both scale together.

**Headline claim:** `maxgap(recalibrated focus map) <= 2/Nt` predicts pass/fail across
all 9 tested configurations. It turns "how much TTD sharing can this system take?" into
an offline geometric check — no channel, no noise, no rate simulation.

### Open — `map_focus_shift.m`
Shared TTD makes the delay layer an array spaced `P*lambda/2`, so it has **grating
lobes every `1/P`**; off `fc` the tilted sub-array pattern may select the wrong lobe,
displacing the focus by an integer multiple of `1/P` — holes, not blur. Signature:
`(realized - designed)*P` clusters on integers. If it does, the fix is lobe
disambiguation and would cost **no extra pilots**. Full reasoning in `THEORY.md` §16.

---

## DISPLACEMENT IS NOT THE FAILURE MODE (`map_focus_shift.m`)

At P=8 the **passing** configuration is scrambled as badly as the failing one:

| P | K | rate | stayed put | mean `|dev|` | maxgap/BW |
|---|---|---|---|---|---|
| 8 | 3 | **42%** | 14.1% | 0.617 | **50.1** |
| 8 | 9 | **100%** | 15.7% | 0.578 | **0.3** |

Only the gap separates them. With the recalibrated lookup, **displacement is
harmless because it is known** — which is exactly why simple recalibration beat a
full match filter in §12: the information was never lost, only mislabeled.

Grating-lobe selection is **not confirmed** (`|k-round(k)|` = 0.09 at P=4 but
0.16-0.24 at P=8, and `|k|=1` is ~0% while `|k|>=2` dominates) and is **moot** —
the statistic is the same for pass and fail, so it is not the discriminator.

### The law now decomposes
    (C1) coverage, DD-BS's own:  K * s >= K0 = 3
    (C2) sharing, new:           P * |theta_t0| * s <= Theta ~ 85
    =>  K_min = K0 * P * |theta_t0| / Theta = 1.12 * P     (was fitted, now derived)

`kspace_map.m` scans K and s independently at P=8 to test it. Both
single-constraint corners (K=9,s=1 and K=3,s=1/3) must FAIL. The resulting heat
map with both boundaries drawn on it is the paper's key figure.

---

## CORRECTION — `P*|theta_t| <= 85` was `K` in disguise (`kspace_map.m`)

Scanning K and `s` **independently** at P=8 kills both constraints of THEORY §17:

| K \ s | 1.00 | 0.50 | 0.333 | 0.250 |
|---|---|---|---|---|
| 3 | 43% | 46% | 39% | 42% |
| 6 | 72% | 87% | 86% | 80% |
| 9 | **100%** | **102%** | 100% | 96% |
| 12 | 101% | 101% | 99% | 100% |

`K=9, s=1` has `P*|theta_t| = 253.8` — 3x over the supposed limit — and gives 100%.
Rows are flat; **only K matters.** The false invariant arose because
`redesign_for_shared` walked the diagonal `s = K0/K`, where `theta_t` is a function
of K. *Never fit an invariant on a locus where the candidate variables are tied.*

**What survives is simpler:** `K_min` depends on **P alone** (3, 5, 9, 18 for
P = 2, 4, 8, 16), independent of `s`, `theta_t` and `B`. So the two knobs separate:
pay `K = ceil(1.12 P)` for the TTD budget (unavoidable), then set `s = K0/K` for
**free** — rate is flat along the row — and collect the 4x delay-range reduction
as a bonus rather than a purchase.

**The maxgap criterion is now 25/25** across every configuration tested, having
outlived four mechanism hypotheses. It is the centerpiece.

### Next — `pilot_spacing_map.m`
The baseline hard-wires the pilot comb to `Delta = 2/K`. Untying K from `Delta`
tests the one hypothesis that explains the null results as well as the positive
ones: `Delta <= 2/P` (comb finer than the sub-array beam) and `K*Delta >= 2` (span).
If `K=8, Delta=0.25` passes while `K=16, Delta=0.5` fails, the fix is **place the
pilots better**, not **use more of them**.

---

## CLOSED — one master curve; non-uniform placement refuted

**Gap statistic fixed, criterion restored 20/20** on the re-run: the two
counterexamples now read 128.7 beamwidths instead of 0.6/0.3. (The 9+16 earlier
configurations were scored with the buggy version; conclusions unaffected, but the
"25/25" figure is withdrawn until they are re-run.)

**The master curve.** Among all configurations satisfying `K*Delta >= 2`, rate
depends on `Delta*P/2` — comb spacing over sub-array beamwidth — and on nothing else:

| `Delta*P/2` | 0.50 | 0.67 | 1.00 | 1.14 | 1.33 | 2.00 |
|---|---|---|---|---|---|---|
| rate | 102-103% | 101% | 100-103% | 98% | 73% | 52-59% |

At 2.00 the rate is 52-59% for K = 4, 6, 8, 12 **and** 16 — a four-fold change in
pilot count moves nothing. **K drops out.** One curve explains the whole table;
this is the paper's key figure.

**`K_min = 7` at P=8** (K=6 → 73%, K=7 → 98%), so `K_min ~ 0.875 P` and the
exchange law is `K * N_TTD ~ Nt`. Earlier values (1.12P, then P) were ladder
artifacts; this one comes from a spacing sweep, not a bracket.

**Non-uniform placement is refuted** (`nonuniform_comb.m`): greedy loses 21-24
points exactly where it matters (K=8: 79% vs 100%), because its offsets leave a
0.624 gap against the 0.25 limit. And it cannot be repaired — the required
resolution `2/P` is set by the sub-array beamwidth, which does not vary with
position, so the optimal fixed-budget comb *is* the uniform one. **The pilots must
be paid; there is no placement trick.**
