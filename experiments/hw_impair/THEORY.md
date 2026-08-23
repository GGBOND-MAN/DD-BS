# TD resolution requirement for DDBS beam training — derivation

> ### ERRATUM (supersedes §7 as first written)
> The Kronecker structure of §7 was originally claimed to cut the **number of TTD
> elements** from `Nt` to `Nt/P + P` (an "8x reduction"). **That claim is wrong and
> is retracted.** `Nt/P + P` counts *distinct delay values*, not physical elements.
> Walking the signal flow shows why: in AoSA sharing the delay is applied **before**
> the fan-out, so `P` antennas genuinely share one element (`G = Nt/P` elements). In
> the Kronecker structure antenna `n = (g,q)` needs `A(n_g) + C(n_q)`; `A` can be
> applied before the fan-out, but each of the `P` antennas then needs its **own**
> physical `C` element — and the `C` values, though numerically identical across
> sub-arrays, act on different signals and cannot share hardware. Element count is
> `G + Nt`, i.e. **more** than a fully-connected array, in either ordering.
> What the structure does reduce is the number of **distinct delay values**
> (`Nt -> Nt/P + P`), which matters for a fixed-line-plus-switch implementation
> (fewer line lengths to fabricate and calibrate) but is a far smaller benefit than
> claimed. The measured gains in §7-§10 are unchanged; only the **cost axis they
> are plotted against** was mislabelled. §7 is rewritten below accordingly.

Analytic backing for the measured naive-vs-hybrid gap. All formulas are validated
against the end-to-end simulation at the bottom.

## 1. Quantization error → array gain

The per-antenna phase of a DDBS beam is `psi_{n,m} = 2*pi*f_m*tau_n + phi_n`, with

    tau_n = (n*d*theta_t - n^2*d^2*alpha_t)/c .

Quantizing the delay with LSB `D` gives `delta_n ~ U[-D/2, D/2]`, hence a phase
error `eps_{n,m} = 2*pi*f_eff*delta_n ~ U[-pi*f_eff*D, pi*f_eff*D]`, where `f_eff`
is the frequency that multiplies the delay error (see §2). Averaging over the
aperture, the realized focusing gain is

    G = |E[exp(j*eps)]| = sinc(f_eff * D),        sinc(x) = sin(pi*x)/(pi*x).

For a target gain `G`, write `x_G = sinc^{-1}(G)` (e.g. `x_0.99 = 0.078`,
`x_0.90 = 0.250`). The resolution requirement is

    **D <= x_G / f_eff .**

Note `G` goes negative past the first sinc null — that is the observed *cliff*,
not a gradual roll-off.

## 2. The two architectures differ only in `f_eff`

**Pure-TD (baseline).** The TTD realizes the whole term `2*pi*f_m*tau_n`, so a
delay error is amplified by the full carrier: `f_eff = f_H = fc + B/2`.

**Hybrid TD-PS (proposed).** Split at the carrier,

    2*pi*f_m*tau_n = 2*pi*fc*tau_n + 2*pi*(f_m - fc)*tau_n .

The first term is frequency-independent and is realized *exactly* by the phase
shifter (mod 2*pi); the TTD only carries the residual, so `|f_m - fc| <= B/2` and

    f_eff = B/2 .

## 3. Delay range, and why it cannot simply be shrunk

The linear term dominates (`R_lin = 140 ns` vs `R_quad = 0.72 ns` at the paper's
settings), so the delay range is `R ≈ Nt*d*|theta_t|/c`. Substituting the paper's
own design equations (29)–(30),

    theta_p + 2*p_M = 1.76*gamma*fL*M/(Nt*B),      theta_t = 1 - (fc/fM)*(theta_p + 2*p_M)

(which reproduce the paper's `theta_t = -32.95` exactly) and `d = c/(2*fc)` gives

    **R = 1.76*gamma*fL*M / (2*fM*B)  ≈ 0.71 * (M/B)  = 145 ns**   — *independent of Nt*.

So the large delay range is **intrinsic to DDBS**: it is set by the angular sweep
rate the scheme needs (eq. 28), i.e. by the OFDM symbol scale `M/B`, not by the
array size. Reducing it is not available as a lever — which is why the precision
requirement, not the range, is the thing to attack.

## 4. Required bits

With `bits = log2(R/D)`:

    B_pureTD  >=  log2( R * f_H / x_G )
    B_hybrid  >=  log2( R * B / (2 * x_G) )

and therefore

    **B_pureTD - B_hybrid = log2( 2*f_H / B )**

— independent of `Nt`, `theta_t`, `R`, *and* of the target gain `G`. It depends
only on the **fractional bandwidth**, and grows as the band narrows:

| B (fc = 30 GHz) | 1 GHz | 2 GHz | 5 GHz | 10 GHz |
|---|---|---|---|---|
| bits saved | 5.9 | 5.0 | **3.7** | 2.8 |

At the paper's `B = 5 GHz`: pure-TD needs ~14.2 bits, hybrid ~10.5 bits (target
`G = 0.90`).

## 5. Validation against the end-to-end simulation

Predicted gain `sinc(f_eff*D)` vs measured average rate (ideal DDBS = 4.771):

| B_td | LSB [ps] | pred pure-TD | meas pure-TD | pred hybrid | meas hybrid |
|---|---|---|---|---|---|
| 16 | 2.1 | 0.992 | 100% | 1.000 | 100% |
| 14 | 8.5 | 0.878 | 94% | 0.999 | 100% |
| 13 | 17.1 | 0.564 | 73% | 0.997 | 100% |
| 12 | 34.2 | *null* | **2%** | 0.988 | **100%** |
| 11 | 68.4 | — | 2% | 0.953 | 99% |
| 10 | 136.8 | — | 3% | 0.819 | 96% |
| 9 | 273.5 | — | 3% | 0.390 | 76% |
| 8 | 547.0 | — | 2% | *null* | 41% |

The model places the pure-TD cliff between 13 and 12 bits and the hybrid cliff
between 9 and 8 bits — both exactly where the simulation shows them.

## 6. What this gives the paper

- A closed-form TD-resolution requirement for DDBS, in the paper's own notation.
- The delay range is intrinsic (`~0.71*M/B`, independent of `Nt`) — so precision,
  not range, is the design lever.
- A hybrid TD-PS split that relaxes the requirement by `log2(2*f_H/B)` bits
  (3.7 here, ~6 at 1 GHz), verified end-to-end, with the phase shifter *not*
  becoming the new bottleneck (3-bit PS still gives ~97%).
- Orthogonal to the OJ-COMS AoSA result, which reduces the TD *count* assuming
  ideal delay values. This work addresses the *precision* dimension that they (and
  all prior DDBS work) assume away, so the two **compose**: the reparameterization
  can be applied on top of their architecture. That composition — not a claim of
  beating them — is the defensible framing.

---

## 7. A second, weaker lever: reducing the number of **distinct delay values**

### Why naive sub-array sharing is expensive for DDBS
Sharing one TTD across `P` antennas leaves a residual `2*pi*(f_m-fc)*(tau_g - tau_n)`.
The intra-group delay spread is `P*d*|theta_t|/c`, and `|theta_t| ~ 33` for DDBS,
so even `P = 2` gives ~1.1 ns of spread — several radians at `B/2`. Measured
focusing gain: `P = 2 -> 0.668`, `P = 4 -> 0.402`, `P = 8 -> 0.238`.

### A Kronecker factorization of the delay profile
**Scope (see the erratum above): this reduces distinct delay *values*, not the
number of physical delay elements.** Split the antenna index exactly as
`n = n_g + n_q` (`n_g` the sub-array offset,
`n_q` the position within it). Then

    tau_n = (n*d*theta_t - n^2*d^2*alpha_t)/c
          = [A(n_g) + C(n_q)]  -  2*n_g*n_q*d^2*alpha_t/c ,

i.e. the **linear term is exactly separable** and only the quadratic **cross term**
resists. Since the DDBS delay is linear-dominated (140 ns vs 0.72 ns here — the
same large `theta_t` that makes the hardware hard is what makes this factorization
near-exact), a two-stage bank of

    Nt/P coarse TTDs  +  P fine TTDs (reused by every sub-array)  =  Nt/P + P

reproduces the delay profile using only `Nt/P + P` **distinct delay values**,
minimised at `P = sqrt(Nt)` -> `2*sqrt(Nt)` values (32 distinct values instead of
256 at `Nt = 256`). Physical element count is **not** reduced — see the erratum.
The practical benefit is therefore limited to fabrication/calibration of a
fixed-delay-line bank, and this should be presented as a secondary observation,
not as a hardware-cost contribution.

### Design rule for P
The residual is bounded by `d_tau <= Nt*P*d^2*|alpha_t| / (2*c)`, so the worst-case
phase error is

    **eps_max = pi * B * Nt * P * d^2 * |alpha_t| / (2*c)**

— it grows **linearly in P**, while the TTD count falls as `Nt/P + P`. Choose
`P = min( sqrt(Nt), 2*c*eps_budget/(pi*B*Nt*d^2*|alpha_t|) )`. At the paper's
settings both terms land near `P = 16` (`eps_max = 1.4 rad`, gain 0.941).

### Measured surface (ideal TD, hybrid PS; ideal full-TTD = 0.975)

| distinct delay values | 128 | 64 | 32 |
|---|---|---|---|
| generic sharing | 0.668 | 0.402 | 0.238 |
| **Kronecker** | **0.975** | **0.973** | **0.941** |

~4x higher gain at equal TD count. And the two contributions **compose**:
Kronecker at 32 distinct values still gives **0.938 at 12-bit**. Note this is a
*distinct-value* axis, not an element-count axis (erratum) — the load-bearing
result on this plot is the **precision** relaxation, which holds for every
architecture.

`P = 16 = sqrt(Nt)` is confirmed optimal: `P = 32` needs *more* TTDs (40) and
drops to 0.853.

## 8. The two contributions are orthogonal — a natural ablation

The corrected sweep shows the `shared` arm is almost **bit-insensitive**
(0.668 at infinite resolution, 0.657 at 11 bits). Once the PS carries the exact
`2*pi*fc*tau_n`, *both* the sharing error and the quantization error scale with
`(f_m - fc)`. So:

* **Contribution A (hybrid TD-PS) is architecture-independent** — it repairs the
  precision axis for `full`, `shared` and `kron` alike.
* **Contribution B (Kronecker) acts on the distinct-value axis** — it is what
  keeps a 32-distinct-value delay profile usable. Per the erratum this is a much
  weaker hardware claim than originally stated, and should be reported as a
  secondary observation.

| architecture | 32 values, 12-bit | axis addressed |
|---|---|---|
| AoSA form, unchanged DDBS parameters | 0.060 | none |
| + hybrid TD-PS (A) | 0.237 | precision |
| + Kronecker (B) | **0.938** | precision + distinct values |

`hw_probe_ablation_rate.m` reproduces this table in the rate domain.

## 9. Positioning — two caveats that must survive into the paper

**Do not claim to beat OJ-COMS.** The `ojcoms` arm applies their *architecture*
(their eqs. (13)-(14): TTD on the sub-array-centre grid, plain per-antenna DDBS
phase shifter that does **not** compensate the intra-group delay) to the
**original** DDBS parameters. Their complete scheme additionally redesigns those
parameters and adds sector-shift, distance interleaving and a larger pilot budget
(~12 vs 3), under which they report near-baseline performance at half the TTD
count. The defensible comparison is against the `shared` arm — their TTD sharing
combined with contribution A, i.e. the strongest fair sharing baseline.

**Prior art on the precision axis concluded the opposite — and that is the point.**
P.-H. Chang and T.-D. Chiueh, "Hybrid beamforming for wideband terahertz massive
MIMO communications with low-resolution phase shifters and true-time-delay,"
*IEEE TWC*, 23(7):8000-8012, 2024 **does** quantize the TTD (its eq. (32), via a
soft-quantization/QAT gradient method) and reports that "*the spectral efficiency
reduction is quite limited*" as the delay interval coarsens.

| | Chang & Chiueh 2024 | this work |
|---|---|---|
| task | hybrid **beamforming** (data precoding) | DDBS **beam training** |
| channel | 3D geometric S-V, AoD/elevation — **far field** | **near field**, spherical, with range |
| delay range | set by the array aperture, **~2 ns** | set by `theta_t ~ -33`, **~140 ns (67x)** |
| TTD LSB studied | `1/(8*fc)` ~ 4.2 ps — already fine | swept to the failure point |
| conclusion | coarse TTD is tolerable | coarse TTD is **catastrophic** |

The two are consistent under the mechanism of §1-§4: at equal *absolute* LSB the
DDBS training waveform needs `log2(67) ~ 6` more bits. Prior studies concluded
"coarse TTD is fine" because their delay range is aperture-limited; the DDBS
training waveform is the outlier — which is exactly why this bottleneck was
missed. State the contrast explicitly: it distinguishes the contribution *and*
explains the gap in the literature.

## 10. Rate-domain ablation (measured) — `hw_probe_ablation_rate.m`

All arms at a matched **32 distinct delay values** (`ojcoms`/`shared`: P=8 ->
Nt/P, which for these two *is* also the element count since sharing applies the
delay before the fan-out; `kron`: P=16 -> Nt/P+P distinct values, element count
unchanged — see erratum). Training and serving beams share one hardware model (same architecture,
same absolute delay LSB). SNR 15 dB, K=3, single path.
Ideal full-TTD DDBS = **4.753**, Perfect CSI = 5.028 bit/s/Hz.

| arch | Inf | 14 | 13 | 12 | 11 | 10 bits |
|---|---|---|---|---|---|---|
| `ojcoms` form | 0.385 | 0.392 | 0.273 | 0.162 | 0.157 | 0.171 |
| `+ hybrid TD-PS` | 0.495 | 0.478 | 0.595 | 0.753 | 0.646 | 0.622 |
| **`+ Kronecker`** | **4.768** | **4.763** | **4.760** | **4.761** | **4.728** | **4.601** |

As a fraction of ideal: ~3-8% / ~10-16% / **97-100%**. The proposal holds
near-ideal rate at a **32-value delay profile and 12-bit resolution**, where both
sharing arms have collapsed. Report the **12-bit** part as the contribution; the
value-count part is secondary.

### Three cautions when reporting this
1. **The `shared` row is non-monotonic in bits** (0.495, 0.478, 0.595, 0.753,
   0.646, 0.622). That is *not* "12-bit is best" — the arm is in a fully broken
   regime where which spurious peak wins is essentially random, so 100-trial
   Monte-Carlo noise dominates. Report it only as "~10-16% of ideal, i.e. broken";
   never quote an ordering within the row.
2. **Rate collapses harder than focusing gain** (gain 0.238 -> rate ~10%): gain is
   an amplitude, power squares it (~-12 dB), and a mis-located argmax additionally
   mis-points the serving beam. Expected, not a bug.
3. **Kronecker's rate does not show its 3.5% focusing-gain loss** (0.941 vs 0.975).
   Legitimate: training is an argmax, so a moderate gain loss still selects the
   right grid point, and the rate is then set by the serving beam, whose 2.1 ns
   delay range the Kronecker structure reproduces almost exactly.

### An explanatory hypothesis (state as such, do not claim it)
This result offers a reason why OJ-COMS needs ~12 pilots: plain sub-array sharing
applied to the *original* DDBS parameters fails, so their scheme must recover it
through parameter redesign, sector-shift and a larger pilot budget. The Kronecker
structure instead keeps near-ideal rate at the original **3** pilots. Since their
full scheme is not reproduced here, this must be offered as an explanation, never
as "we beat them with fewer pilots".

## 11. The reparameterization generalizes: it attenuates **any** delay error

The composition run (`hw_probe_composition.m`, P=2, 128 shared TTDs,
**N_iter = 500**, ideal full-TTD = 4.754) gave:

| architecture | Inf | 15 | 14 | 13 | 12 | 11 | 10 bits |
|---|---|---|---|---|---|---|---|
| AoSA form (their split) | 2.730 | 2.616 | 2.511 | 1.276 | 0.096 | 0.099 | 0.106 |
| **AoSA + reparameterization** | **4.192** | 4.277 | 4.238 | 4.164 | 4.173 | 4.199 | 4.035 |

As a fraction of ideal: their split runs 57% -> 2%; ours 88% -> 85%. At 12 bits
the two differ by **43x**.

Two things follow.

**(a) The precision claim, on a fixed architecture.** Staying within 5% of each
arm's *own* infinite-resolution rate: their split needs **15 bits** (14 already
fails at 2.511 vs a 2.594 threshold), the reparameterized version is still within
5% at the sweep floor of **10 bits** — a saving of **>= 5 bits** on an identical
AoSA array. This is the composition claim, and it needs no reproduction of their
scheme.

**Exact figure (extended sweep to 7 bits).** With the sweep extended, our arm's
own cliff is found and the saving becomes exact:

| architecture | Inf | 15 | 14 | 13 | 12 | 11 | 10 | 9 | 8 | 7 |
|---|---|---|---|---|---|---|---|---|---|---|
| their split | 2.664 | 2.647 | **2.302** | 1.333 | 0.077 | 0.094 | 0.093 | 0.098 | 0.086 | 0.082 |
| + reparam. | 4.185 | 4.329 | 4.191 | 4.285 | 4.299 | **4.072** | 3.921 | 2.787 | 1.158 | 1.010 |

Minimum resolution within 5% of the same arm's infinite-resolution rate:
**15 bits** for their split, **11 bits** with the reparameterization —

    saving = 4 bits,   analytic prediction = log2(2*f_H/B) = log2(13) = 3.70 bits.

An independent cross-check of the theory against simulation: the predicted 3.70
rounds to the measured 4. The cliff locations shift consistently too — theirs
between 13 and 12 bits, ours between 10 and 9.

**Precision of that integer.** The 10/11-bit boundary sits inside the Monte-Carlo
noise: at `N_iter = 100` the 10-bit point reads 3.921 and fails the 5% band
(giving 11), while at `N_iter = 500` it read 4.035 and passes (giving 10). Quote
the result as **~4 bits (4-5)**, consistent with the analytic 3.70, and re-run the
full extended sweep at `N_iter = 500` for the published figure.

**(b) It is not only about quantization.** At *infinite* resolution the
reparameterization already wins, 2.730 -> **4.192 (+54%)**. There is no
quantization error there, so the gain must come from elsewhere — and it does. The
error model of §1-§2 never referred to quantization specifically:

    any delay error d_tau  ->  phase error  2*pi * f_eff * d_tau

`d_tau` can come from quantization, from **sub-array sharing** (`tau_g - tau_n`),
or from fabrication jitter. The reparameterization lowers `f_eff` from `f_H` to
`B/2` for **all of them alike**, a ~13x attenuation at these parameters. So it
also makes AoSA sharing itself substantially more viable — a stronger and more
general statement than "it fixes quantization", and the one the paper should make:

> On a TTD-PS array, moving the frequency-independent term `2*pi*fc*tau_n` into
> the phase shifter attenuates the beam-degradation caused by *any* delay error —
> quantization, sub-array sharing, or fabrication — by the factor `2*f_H/B`.

**Residual gap.** The reparameterized arm reaches 4.19 of the ideal 4.754 (88%);
the remaining 13% is the sharing error that survives the `B/2` scaling. Honest to
state: the reparameterization *attenuates* the sharing error, it does not remove it.

**Monte-Carlo caution.** At `N_iter = 100` the flat arm fluctuated by several
percent (a 12-bit point read 4.362, above its own infinite-resolution 4.132). At
`N_iter = 500` the spread tightens to ~3% and that artefact is gone, but a mild
wiggle remains (15 bits reads 4.277 vs 4.192 at infinite resolution). The arm is
flat within noise — do not read an ordering into those wiggles.

---

## 13. PRIOR-ART FINDING #2 — the resolution fix is also prior art, but a bigger DDBS-specific gap opens

The paywalled paper has been obtained and read:

> A. Najjar, M. El-Absi, T. Kaiser, "Hybrid Delay-Phase Precoding in Wideband
> UM-MIMO Systems Under True Time Delay and Phase Shifter Hardware Limitations,"
> *IEEE TWC* 23(7):7246-..., 2024.

Its **Sec. V-B, "True Time Delay Resolution Constraint"**, analyses exactly the
error "*that occurs due to the adoption of finite-resolution TTDs ... from rounding
the required time delay value to an integer multiple of the time delay resolution
`t_res`*", introduces an auxiliary variable via `phi_k^n = 2*pi*rho_k^n*tau_k^n`,
minimises the squared phase error over all subcarriers, and obtains

> "As a result of (37), **`rho_k^n = fc`**, and, accordingly, the corresponding
> **equalization phase is `phi_k^n = 2*pi*fc*tau_k^n`**. The equalization phase
> `phi_k^n` is added to the corresponding subarray..."

That is the proposed mechanism, derived for the proposed motivation (finite TTD
resolution), by a least-squares argument that is cleaner than the heuristic used
here. **Contribution A is therefore comprehensively prior art**: the split itself
in Dai/Tan TWC 2022, and its use against quantization error in Najjar et al. 2024.
It must be cited and cannot be claimed in any form.

### But the same reading exposes a larger, unreported, DDBS-specific problem

Delay *ranges* actually realizable, as surveyed by that paper and by Dai/Tan:

| source | max delay range |
|---|---|
| THz SiGe/mHEMT TTD circuits [28],[30],[31] of Najjar et al. | **1.47 - 6.64 ps** |
| mmWave TTDs cited by Dai/Tan | **400 - 508 ps** |
| Najjar et al. simulation settings | `t_max` = 10 / 28 / 340 ps |
| aperture-limited near-field focusing beam (Nt=256, fc=30 GHz) | ~4.25 ns |
| **DDBS training waveform (this work)** | **140 ns** |

DDBS needs **~280x the best realizable device** and **~33x even the
aperture-limited focusing delay** — about **21 m of transmission line per
antenna**, for 256 antennas. **The baseline never states this number.** It says
only that "*only K (usually a small number such as 2,3) time delay values need to
be realized*" and proposes fixed microstrip lines or waveguides, without computing
their required length.

And per Sec. 3 the range is `R = 1.76*gamma*fL*M/(2*fM*B) ~ 0.71*(M/B)`,
**independent of `Nt`** — so unlike every other delay requirement in the
literature, it does **not** shrink with a smaller array.

### The tradeoff this opens (verified numerically)

The range is set by the angular sweep rate (baseline eqs. (28)-(30)). Reducing the
sweep rate by `S` (so each pilot lays down fewer strips and ~`S` times more pilots
are needed for the same coverage) scales the intercept and hence the range:

| sweep reduction `S` | pilots | `theta_t` | delay range | vs a 508 ps device |
|---|---|---|---|---|
| 1 (baseline) | 3 | -32.95 | **140.6 ns** | 277x over |
| 2 | ~6 | -15.98 | 68.2 ns | 134x |
| 5 | ~15 | -5.79 | 24.7 ns | 49x |
| 10 | ~30 | -2.40 | 10.2 ns | 20x |
| 20 | ~60 | -0.70 | **3.0 ns** | ~6x |

`R ∝ 1/S`, flooring at the aperture-limited ~4 ns once `theta_t -> 0`.

**So DDBS's headline "3 pilots" is paid for in TTD delay range, at a rate nobody
has quantified.** That is a genuinely open, DDBS-specific, and consequential
contribution — and it is *orthogonal* to the resolution question that turned out
to be prior art:

1. DDBS's required delay range is 140 ns, never stated in the baseline;
2. it scales as `0.71*M/B`, **independent of `Nt`** — a qualitatively different law
   from every aperture-limited delay requirement in the DPP literature;
3. it exceeds realizable TTD hardware by 2-3 orders of magnitude;
4. it trades against pilot overhead as `R ∝ 1/S`, with a floor at the aperture
   limit — giving designers the first realizable operating point for DDBS.

---

## 14. MEASURED — parameter redesign makes shared-TTD DD-BS work, and the tradeoff obeys a clean law

`redesign_for_shared.m`, Nt=256, fc=30 GHz, B=5 GHz, M=1024, SNR=15 dB, N_iter=60,
single-path LoS, decoded with the **recalibrated lookup** (`actual_focus`):

| `s` | `K` | `theta_t` | delay range | 256 TTD | 128 TTD | 64 TTD | 32 TTD |
|---|---|---|---|---|---|---|---|
| 1.00 (baseline) | 3 | 31.73 | 134.8 ns | 4.740 | 4.689 | 3.665 | **1.683** |
| 0.50 | 6 | 15.86 | 67.4 ns | 4.788 | 4.771 | 4.686 | 4.127 |
| 0.25 | 12 | 7.93 | 33.7 ns | 4.744 | 4.733 | 4.721 | **4.767** |

### What the table says

1. **Read across the top row** — the baseline design point collapses under sharing:
   64 TTD → 77%, 32 TTD → **35%** of ideal. This is the negative answer to
   research question 1 (*does the original architecture/algorithm survive limited
   TTD?* — **no**), already established in §12 and reconfirmed here.
2. **Read the last column downward** — at a *fixed* 32-TTD budget, redesigning the
   sweep rate takes the rate from 1.683 → 4.127 → **4.767**, i.e. from 35% to
   **100.6%** of the full-per-antenna ideal (4.740). The loss is not intrinsic to
   sharing; it is intrinsic to *using the baseline's parameters* under sharing.
3. **Read the first column downward** — 4.740 / 4.788 / 4.744, flat. Slowing the
   sweep costs **no rate at all** when TTDs are plentiful. So `s` is a pure
   *hardware-for-pilots* exchange, not a rate sacrifice — the comparison is clean.
4. **The price is pilots, and it is bounded.** K goes 3 → 12. In the baseline's own
   Table I the competing schemes need 10 (near-field rainbow), 140 (hierarchical),
   266 (two-stage), 2560 (exhaustive). K=12 is at parity with near-field rainbow
   and still 12x-200x below the rest — while using **8x fewer TTD elements** than
   any of them.

### Two hardware savings at once

The redesign also shrinks the delay **range** (§13's unrealizability problem):
134.8 ns → 33.7 ns, a 4x reduction, because the range is set by `|theta_t|` and
`|theta_t| ∝ s`. Total delay-line budget (elements x range):

    baseline   256 x 134.8 ns = 34509 ns-elements
    redesigned  32 x  33.7 ns =  1078 ns-elements     -> 32x less

### The design law: `K_min ~ 1.5 P`, i.e. `K * N_TTD ~ 1.5 Nt`

Taking the smallest K that reaches the ideal rate in each column:

| N_TTD | 256 | 128 | 64 | 32 |
|---|---|---|---|---|
| K_min | 3 | 3 | 6 | 12 |
| K/P | 3 | 1.5 | 1.5 | 1.5 |

`K/P >= 1.5` holds at **every** entry of the table, passing and failing alike
(P=4,K=3 → K/P=0.75 fails at 77%; P=8,K=6 → 0.75 fails at 87%). So the frontier is
the hyperbola `K * N_TTD = 1.5 * Nt = 384`.

**Where the constant comes from.** Antennas in a group of P share one delay, so the
residual phase is `2*pi*(f_m - fc)*(tau_g - tau_n)` with `|f_m - fc| <= B/2`, and the
intra-group delay spread is dominated by the linear DDBS term, `d_tau ~ P*d*|theta_t|/c`.
Holding the worst-case residual phase below a fixed `phi_max`, with `d = c/(2 fc)`:

    P * |theta_t|  <=  4 * fc * phi_max / (pi * B)                      (INVARIANT)

and since the redesign ties `|theta_t| = |theta_t0| * K0 / K`,

    K / P  >=  (pi * B / (4 * fc * phi_max)) * |theta_t0| * K0   ∝  B / fc.

The measured `K/P = 1.5` at `B/fc = 1/6` pins `phi_max ~ 8.6 rad` and predicts
`K/P = 9 * (B/fc)`. **This prediction is falsifiable and untested**: at P=8, K_min
should be 6 / 12 / 24 for B = 2.5 / 5 / 10 GHz. `frontier_ttd_pilot.m` runs exactly
that check (plus a finer frontier including P=16). Until it passes, the constant is
empirical at one operating point and must be reported as such.

### Status against the user's three research questions

| Q | status |
|---|---|
| 1. Does the original architecture/algorithm survive limited TTD? | **Answered: no.** 32 TTD → 35% of ideal. |
| 2. Design something that approaches ideal under limited TTD | **Answered, two independent levers, both measured**: (i) recalibrated lookup — free, algorithm-only, 32 TTD 10% → 38%; (ii) parameter redesign — 32 TTD → 100.6% at K=12. Together they are the paper. |
| 3. Abandon DD-BS for a new architecture | not started; only worth doing if (2) turns out to be beatable. |

### What is genuinely new here vs. §13's prior art

§13 killed contribution A (the TD/PS split) as prior art. **This is a different
claim** and none of the surveyed work states it:
- the shared-TTD failure of DD-BS is not a resolution problem but a **parameter-choice**
  problem (Najjar et al. treat resolution; OJ-COMS changes the architecture but keeps
  ideal per-group delays and does not report a K-vs-TTD frontier);
- the exchange rate `K * N_TTD = const * Nt` with `const ∝ B/fc`;
- the resulting joint reduction of **TTD count and TTD delay range** — the latter
  being the §13 unrealizability problem, which the same knob attacks.

---

## 15. MEASURED + FALSIFIED — the frontier is `K_min = ceil(1.12 P)`, and it does **not** scale with `B`

`frontier_ttd_pilot.m`, same setup as §14, N_iter=30, pass = 95% of that
configuration's own full-TTD ideal.

### PART 1 — the frontier, at finer K resolution

| P | N_TTD | K_min | K/P | `theta_t` | range | `P*|theta_t|` |
|---|---|---|---|---|---|---|
| 2 | 128 | 3 (floor) | 1.50 | 31.73 | 134.8 ns | 63.5 |
| 4 | 64 | 5 | 1.25 | 19.04 | 80.9 ns | 76.1 |
| 8 | 32 | 9 | 1.12 | 10.58 | 44.9 ns | 84.6 |
| 16 | **16** | 18 | 1.12 | 5.29 | 22.5 ns | 84.6 |

Failing points bracket each: P=4,K=3 → 74%; P=8,K=3 → 42%, K=6 → 89%;
P=16,K=6 → 58%, K=12 → 80%. P=2's `K_min` is limited by the `K >= K0 = 3` floor,
not by sharing.

**`K/P` is not constant — `P*|theta_t|` is.** The invariant is

    P * |theta_t|  <=  Theta ~ 85          (dimensionless)

and since the redesign sets `|theta_t| = |theta_t0| K0 / K`, this is

    K_min = ceil( P * |theta_t0| * K0 / Theta ) = ceil(1.12 * P)

which reproduces **all four** measured points exactly (P=2→3 at the floor, 4→5,
8→9, 16→18). §14's `K ~ 1.5P` was a coarse-ladder artifact; supersede it.

### PART 2 — the `B/fc` prediction is FALSIFIED

At P=8, sweeping B with each point judged against its own full-TTD ideal:

| B | ideal rate | K=3 | K=6 | K=9 | K_min | predicted K/P |
|---|---|---|---|---|---|---|
| 2.5 GHz | 3.834 | 42% | 93% | 98% | 9 | 0.75 |
| 5 GHz | 4.742 | 43% | 86% | 100% | 9 | 1.50 |
| 10 GHz | 4.449 | 37% | 94% | 106% | 9 | 3.00 |

`K_min = 9` at every bandwidth — a 4x change in B moves it not at all, against a
predicted 4x change. **§14's derivation is wrong.** The worst-case residual phase
`pi B P d |theta_t| / c` cannot be the binding constraint, because it is
proportional to B. (Sanity check on its magnitude: at the P=8 frontier it equals
11 rad at B=5 GHz and 22 rad at B=10 GHz — both far past any sensible `phi_max`,
yet both pass. The bound is simply not what is being enforced.)

### Replacement mechanism — coverage truncation, in which both `fc` and `B` cancel

1. The intra-group residual `2*pi*(f_m-fc)*(tau_g-tau_n)` is **linear in the
   intra-group index**, because the DDBS delay is linear-dominated. A linear phase
   across a sub-array is a **tilt**, not a loss: each sub-array pattern is steered
   off by `d(sin theta) = -(f_m - fc) theta_t / f_m`.
2. The tilt leaves the sub-array beamwidth (`~2/P`) once
   `|f_m - fc| > eta * c/(P d |theta_t|) = eta * 2 fc /(P |theta_t|)`, so only that
   slice of the band still forms a usable beam (`eta ~ 0.443` at 3 dB).
3. DDBS sweeps at `d(theta_m)/df ~ -theta_p/fc`, so the **usable angular coverage
   per pilot** is

       2 * (theta_p/fc) * eta * 2 fc /(P |theta_t|)  =  4 eta theta_p / (P |theta_t|)

   — `fc` cancels, and **`B` never enters**. That is exactly the measurement.
4. Covering `Theta_req` with K pilots therefore needs

       K  >=  Theta_req * P * |theta_t| / (4 eta |theta_p|)                     (*)

   With `Theta_req = 2 sin 60 = 1.732`, `|theta_t0| = 31.73`, `eta = 0.443` and
   `|theta_p,eff|` between 30 (`theta_2`) and 36.78 (`theta_p + 2 p_M`):
   `K_min = (0.84 .. 1.04) P` against the measured **1.12 P** — a 8-25% match from
   a first-principles argument with no fitted quantity.

`diag_shared_mechanism.m` tests steps 2-3 directly (measured usable bandwidth vs
`2 fc/(P|theta_t|)`, and the largest coverage hole among surviving subcarriers).

### The consequence — `theta_t` and `theta_p` enter (*) as a RATIO

§14 scaled them **together**, so sharing robustness had to be bought with pilots.
But they do different jobs:

| parameter | role | wants to be |
|---|---|---|
| `theta_p` | sweep rate — strips laid per pilot | **large** (few pilots) |
| `theta_t` | (i) DC offset re-centring the sweep, `~ -theta_p` | **small** (sharing) |
| | (ii) per-pilot increments `2(s-1)/K`, spanning only ~2 | already small |

Only role (i) is large, and by (*) only role (i) causes the failure. Its own
sharing budget is `P * 2 = 16` at P=8 — a factor of 5 inside `Theta ~ 85`.

**So the `K_min = 1.12 P` tradeoff is a property of how DDBS re-centres its sweep,
not of TTD sharing itself.** If the sweep can be re-centred without a large TTD
ramp, `K = 3` and 32 TTDs become simultaneously achievable — the tradeoff
disappears rather than being traded along. `decouple_theta.m` scans `theta_t` and
`theta_p` independently to test this, reporting realized coverage next to rate so
that a beam which merely walked off the user region cannot register as a win.

That is the concrete form of research question 3: **not a new array, but a new way
to place the DDBS sweep** — and it is reached from measurement, not assumption.

---

## 16. MEASURED — the mechanism is COVERAGE HOLES, and it yields a deterministic design test

Two probes, two verdicts. One hypothesis died, one criterion survived and is
stronger than the closed form it replaces.

### 16.1 `diag_shared_mechanism.m` — M1 refuted, M2 confirmed

| P | K | `P*|theta_t|` | BW_meas | BW_pred (M1) | maxgap / (2/Nt) | rate |
|---|---|---|---|---|---|---|
| 1 | 3 | 31.7 | 4.55 GHz | 5.00 GHz | 0.9 | ideal |
| 2 | 3 | 63.5 | 4.55 GHz | 0.84 GHz | 0.9 | 99% |
| 4 | 3 | 126.9 | 3.29 GHz | 0.42 GHz | **20.2** | 74% |
| 4 | 6 | 63.5 | 4.31 GHz | 0.84 GHz | 0.3 | 100% |
| 8 | 3 | 253.8 | 2.64 GHz | 0.21 GHz | **53.6** | 42% |
| 8 | 6 | 126.9 | 4.14 GHz | 0.42 GHz | **10.1** | 89% |
| 8 | 9 | 84.6 | 5.00 GHz | 0.63 GHz | 0.3 | 100% |
| 8 | 12 | 63.5 | 4.67 GHz | 0.84 GHz | 0.3 | 100% |
| 16 | 6 | 253.8 | 4.49 GHz | 0.21 GHz | **25.1** | 58% |
| 16 | 12 | 126.9 | 4.64 GHz | 0.42 GHz | **5.2** | 80% |
| 16 | 18 | 84.6 | 5.00 GHz | 0.63 GHz | 0.3 | 99% |

**(M1) gain truncation is REFUTED.** Measured usable bandwidth is **5-25x wider**
than the sub-array-factor prediction (2.64 GHz vs 0.21 GHz at P=8,K=3). The
sub-array tilt does not remove the band; §15's "usable bandwidth" step is wrong,
and with it the `4 eta theta_p /(P |theta_t|)` coverage formula and (*).

**(M2) coverage holes are CONFIRMED, with perfect discrimination.** The largest
interior gap of the recalibrated focus map, measured in array beamwidths `2/Nt`:

    maxgap <= 1 beamwidth  ->  99-100% of ideal   (4 of 4 configurations)
    maxgap >= 5 beamwidths ->  42-89% of ideal    (5 of 5 configurations)

No overlap, and the transition is a cliff (0.3 → 5.2), not a slope. The `span`
statistic, by contrast, reads 1.04 in **every** row including total failure — a
few stray foci still reach the edges. **Span is useless; the gap is everything.**
(`decouple_theta.m` has been corrected to report `maxgap/(2/Nt)` instead of span.)

### 16.2 `decouple_theta.m` — the ratio design equation is refuted

P=8 (32 TTDs), K=3, scanning `theta_t` and `theta_p` independently. Reference:
full TTD 4.742; P=8 with baseline params 2.042 (43%).

| `s_t` | `s_p` | `P*|theta_t|` | rate | % ideal |
|---|---|---|---|---|
| 1 | 1 | 253.8 | 2.093 | 44% |
| 0.5 | 0.5 | 126.9 | **2.186** | **46%** (best) |
| 0.25 | 1 | 63.5 | 1.416 | 30% |
| 0.125 | 1 | 31.7 | 1.382 | 29% |
| 0.0625 | 1 | 15.9 | 0.576 | 12% |

**Nothing wins.** The ceiling across the whole scan is 46%, and shrinking
`theta_t` at fixed `theta_p` makes it monotonically **worse** — the opposite of
what (*) predicts. Two consequences:

1. **`P*|theta_t| <= 85` is necessary but NOT sufficient.** At `s_t = 0.25` the
   invariant is satisfied (63.5) and the rate is 30%. It only holds along the
   `redesign_for_shared` locus where `theta_t` and `theta_p` scale **together**.
2. **Research question 3 in the "just re-centre the sweep" form is closed.** The
   `theta_t` offset is not a removable nuisance; it is load-bearing. Shrinking it
   alone destroys the beam family rather than freeing it.

### 16.3 What actually survives — and it is a better contribution than the closed form

Both mechanism hypotheses (§15 phase bound, §16 ratio equation) are dead. What is
established by measurement, and is not weakened by their death:

| # | claim | evidence |
|---|---|---|
| 1 | DD-BS collapses under shared TTD at its own parameters | 32 TTD, K=3 → 42-43% of ideal |
| 2 | Recalibrated lookup recovers part of it for free (algorithm only) | §12: 32 TTD 10% → 38% |
| 3 | Joint parameter redesign recovers ~100% | §14: 32 TTD, K=9-12 → 100% |
| 4 | The frontier is `K_min = ceil(1.12 P)`, i.e. `P*|theta_t| <= ~85` | 4/4 points exact |
| 5 | The frontier is **bandwidth-independent** | `K_min = 9` at B = 2.5/5/10 GHz |
| 6 | **A deterministic design test**: `maxgap(recalibrated focus map) <= 2/Nt` | 9/9 configurations, no overlap |

Claim 6 is the one to lead with. It converts "how much TTD sharing can this
system take?" from a Monte-Carlo rate question into a **closed geometric check on
an offline-computable table** — no channel model, no noise realization, no rate
simulation. That is a stronger and more reusable statement than a closed-form
threshold would have been, and it is what makes claims 4-5 usable by a designer.

### 16.4 The one structural hypothesis still standing — `map_focus_shift.m`

Sharing makes the TTD layer an array of `Nt/P` elements spaced `P*d = P*lambda/2`,
so its array factor has **grating lobes every `1/P` in `sin(theta)`**. The
per-antenna PS forms the sub-array pattern that selects the right lobe; off `fc`
that pattern tilts and a **different lobe wins**, displacing the realized focus by
an integer multiple of `1/P`. This produces *holes* rather than blur — which is
what was measured, and which neither dead hypothesis explains.

Falsifiable signature: `(realized - designed) * P` must cluster on **integers**.
`map_focus_shift.m` measures exactly that, plus which designed foci fall into the
biggest hole and where they end up instead.

- clusters on integers with `|k| >= 1` → grating-lobe selection; the fix is lobe
  disambiguation at the sub-array level, and it would be **pilot-free** — a real
  answer to research question 2 that costs no overhead;
- spread uniformly (`|k - round(k)| ~ 0.25`) → the map is smoothly warped, the
  structural route closes, and claims 1-6 stand as the paper.

---

## 17. MEASURED — displacement is NOT the failure mode, and the law splits into two constraints

`map_focus_shift.m`. The decisive comparison is between the two P=8 rows, one of
which passes and one of which fails:

| P | K | rate | stayed put | `|k|>2` | mean `|dev|` | maxgap/BW |
|---|---|---|---|---|---|---|
| 1 | 3 | ideal | 98.3% | 0.0% | 0.036 | 0.6 |
| 4 | 3 | 74% | 28.5% | 31.7% | 0.578 | 20.2 |
| 8 | 3 | **42%** | 14.1% | 61.1% | 0.617 | **50.1** |
| 8 | 6 | 89% | 10.0% | 72.2% | 0.709 | 10.1 |
| 8 | 9 | **100%** | 15.7% | 56.7% | 0.578 | **0.3** |

### 17.1 The finding

**The passing configuration is scrambled just as violently as the failing one.**
At P=8, K=9 — which delivers 100% of ideal — only 15.7% of foci stay put, 56.7%
move by more than two grating steps, and the mean displacement is 0.578 in
`sin(theta)`, all statistically indistinguishable from the K=3 case that delivers
42%. The only column that separates them is the gap: **0.3 vs 50.1 beamwidths**.

So "sharing displaces the beams" is **not** the failure mode. With the
recalibrated lookup, displacement is *harmless because it is known* — the beam may
point anywhere as long as the table says where. What fails is **tiling**.

This retro-validates §12's recalibration result and explains it: the baseline's
decision rule assumes designed focus locations, and under sharing 85% of them are
wrong by more than a grating step. That is why simple recalibration beats a full
match filter — the information was never lost, only mislabeled.

### 17.2 Grating-lobe selection — not confirmed, and moot

`|k - round(k)|` (0 = on a lobe, 0.25 = no structure): 0.093 at P=4, but 0.156 /
0.243 / 0.185 at P=8 — partial structure at P=4, essentially none at P=8. And
`|k| = 1` is ~0% everywhere while `|k| >= 2` dominates, which is the opposite of
what adjacent-lobe selection predicts. Weakly supported at best. **More
importantly it is moot**: the statistic is the same for passing and failing
configurations, so whatever produces the displacement is not the discriminator.
Hypothesis retired; the sub-array-disambiguation route (a pilot-free fix) closes
with it.

### 17.3 What this buys — the law now DECOMPOSES

The frontier `K_min = 1.12 P` was fitted. It is now the intersection of two
independent constraints:

    (C1) COVERAGE   K * s  >=  K0 = 3
         DD-BS's OWN constraint, not new: the sweep per pilot is proportional to
         theta_p ~ s, so total designed coverage is K*s and must be preserved.

    (C2) SHARING    P * |theta_t0| * s  <=  Theta ~ 85
         The new one: above this the displaced foci stop tiling and holes open.

Eliminating `s`:

    K_min = K0 * P * |theta_t0| / Theta = 3 * P * 31.73 / 85 = 1.12 * P

— the measured frontier, now **derived from two constraints** rather than fitted,
with only `Theta` left empirical. And `Theta` never has to be known in closed
form: the maxgap test of §16 computes the C2 boundary offline for any
`(Nt, P, theta_t)`.

This also explains, in retrospect, why the earlier probes behaved as they did:

| probe | what it did | why it gave what it gave |
|---|---|---|
| `redesign_for_shared` | walked the diagonal `K*s = K0` | C1 held by construction; it measured C2 alone |
| `decouple_theta` | broke the tie in the `theta` direction at K=3 | C1 violated → 46% ceiling, no matter what C2 did |
| `frontier_ttd_pilot` | swept P along the diagonal | traced the C1∩C2 corner, i.e. `K_min` |

### 17.4 The falsifiable test — `kspace_map.m`

Scan K and `s` as an **independent 2-D grid** at P=8. The pass region must be
exactly the intersection of the two half-planes. Two single-constraint corners
decide it:

| point | C1 | C2 | prediction |
|---|---|---|---|
| K=9, s=1 | ok (9≥3) | **violated** (253.8 > 85) | FAIL |
| K=3, s=1/3 | **violated** (1 < 3) | ok (84.6 ≤ 85) | FAIL |
| K=9, s=1/3 | ok | ok | PASS |

**Both single-constraint corners must fail.** If either passes, that constraint is
not real and the decomposition is wrong. The resulting (K, s) heat map with the
two boundary curves drawn over it is the paper's key figure.

---

## 18. REFUTED — `P*|theta_t| <= 85` was `K` in disguise. The law depends on `P` alone.

`kspace_map.m`, P=8 (32 TTDs), K and `s` scanned independently. Rate as % of ideal:

| K \ s | 1.00 | 0.50 | 0.333 | 0.250 | maxgap/BW |
|---|---|---|---|---|---|
| 3 | 43% | 46% | 39% | 42% | ~50 |
| 6 | 72% | 87% | 86% | 80% | ~9-10 |
| 9 | **100%** | **102%** | 100% | 96% | 0.3-0.6 |
| 12 | 101% | 101% | 99% | 100% | 0.3-0.6 |

### 18.1 Both constraints of §17 are dead

- **C2 (`P*|theta_t| <= 85`) is REFUTED.** `K=9, s=1` has `P*|theta_t| = 253.8`,
  three times over the supposed limit, and delivers **100%** (maxgap 0.6).
  `K=12, s=1` likewise, at 101%.
- **C1 (`K*s >= K0`) is REFUTED.** `K=9, s=0.25` has `K*s = 2.25` and passes at
  96%; meanwhile `K=3, s=1` satisfies C1 exactly and fails at 43%.

**Rows are flat, columns are everything.** `s` — and therefore `theta_t` — is at
most a second-order effect (its one visible contribution is +15 points in the
transition row, K=6: 72% → 87% as `s` goes 1 → 0.5).

**How the false invariant arose.** `redesign_for_shared` only ever walked the
diagonal `s = K0/K`, where `theta_t = theta_t0 * K0/K` is a *function of K*. So
"`P*|theta_t| = 84.6` at the frontier" and "`K = 9` at the frontier" were the same
statement, and I read the wrong one as causal. §15 and §17's decomposition are
withdrawn. **The lesson is methodological and belongs in the paper's own design:
never fit an invariant on a locus where the candidate variables are tied.**

### 18.2 What survives, and it is simpler than what died

    K_min depends on P alone:  P = 2, 4, 8, 16  ->  K_min = 3, 5, 9, 18
    independent of s (this table), of theta_t (same), and of B (frontier PART 2)

Equivalently **`K * N_TTD ~ 1.12 * Nt`**, the headline exchange law — unchanged in
value, but now known to be a statement about *pilots vs sharing factor*, with the
sweep design playing no part.

**This makes the practical recipe cleaner, not weaker**, because the two knobs
separate completely:

1. TTD budget fixes `P`; pay `K = ceil(1.12 P)` pilots. **Unavoidable.**
2. *Given* that K, set `s = K0/K`. Rate is flat along the row, so this is **free**,
   and it shrinks the delay range 4x (134.8 → 33.7 ns) — §13's unrealizability
   problem, solved at no cost rather than paid for in pilots.

§14's headline ("32 TTDs at 100% of ideal for K=12") stands; only its *causal
story* changes — the pilots buy sharing tolerance, and the range reduction is a
bonus that rides along, not the thing being bought.

### 18.3 The maxgap criterion is now 25/25

Across §16's 9 configurations and these 16, `maxgap <= 1 beamwidth` ⇔ pass, with
**no counterexample in either direction** — including the cases that refuted every
mechanism hypothesis. It has outlived four theories. It is the paper's centerpiece.

### 18.4 The mechanism that fits ALL of it — a pilot comb vs the sub-array beam

The baseline places its pilots by stepping the TTD intercept,
`theta_t,s = theta_t - 2(s-1)/K`, so the pilots form a **comb of spacing
`Delta = 2/K` spanning ~2**. Sharing groups P antennas, whose sub-array pattern has
beamwidth **`~2/P`**. A hole opened by sharing gets filled only if some pilot lands
inside it, i.e. only if the comb is finer than the sub-array beam:

    (R1) RESOLUTION   Delta <= 2/P      ->  with Delta = 2/K:   K >= P
    (R2) SPAN         K * Delta >= 2    ->  the strips must still tile

Neither expression contains `B`, `theta_t` or `s` — which is precisely why `K_min`
was measured to be independent of all three. Together they give `K >= P` against a
measured `1.12 P`. **This is the first hypothesis that explains the null results
as well as the positive one**, which is why it is worth one more test.

### 18.5 The test — `pilot_spacing_map.m`

The baseline hard-wires `Delta = 2/K`. The script generates the pilots one at a
time so **K and `Delta` scan independently**, and the pass region must be the
corner `{Delta <= 2/P}` ∩ `{K*Delta >= 2}`:

| K | `Delta` | span | resolution | prediction |
|---|---|---|---|---|
| 4 | 0.50 | ok | **too coarse** | FAIL |
| 4 | 0.25 | **too narrow** | ok | FAIL |
| **8** | **0.25** | ok | ok | **PASS** — the `K = P` corner |
| 16 | 0.125 | ok | ok | PASS |
| 16 | 0.50 | ok | **too coarse** | FAIL |

**If `K=8, Delta=0.25` passes while `K=16, Delta=0.5` fails, pilot COUNT is not
the variable — pilot RESOLUTION is.** The design fix then becomes *place the
pilots better* rather than *use more of them*, which is strictly cheaper and
would be a second algorithm-side contribution alongside the recalibrated lookup.
If instead only K matters and `Delta` does nothing, the comb model joins the other
four and `K_min ~ 1.12 P` stands as a purely empirical law — still publishable,
with the maxgap criterion as its operational form.

---

## 19. CONFIRMED (20/20) — the variable is pilot RESOLUTION, not pilot count

`pilot_spacing_map.m`, P=8, `s=1` throughout, K and comb spacing `Delta` untied.
Predictions were registered in §18.5 **before** the run.

| K | `Delta` | span `K*Delta` | `Delta*P/2` | rate | % ideal | predicted |
|---|---|---|---|---|---|---|
| 4 | 0.500 | 2.00 | 2.00 | 2.616 | 55% | FAIL ✓ |
| 4 | 0.250 | 1.00 | 1.00 | 2.853 | 60% | FAIL ✓ |
| 4 | 0.125 | 0.50 | 0.50 | 1.956 | 41% | FAIL ✓ |
| 6 | 0.333 | 2.00 | 1.33 | 3.484 | 73% | FAIL ✓ |
| 6 | 0.250 | 1.50 | 1.00 | 3.749 | 79% | FAIL ✓ |
| **8** | **0.250** | 2.00 | 1.00 | **4.728** | **100%** | **PASS ✓** |
| 8 | 0.500 | 4.00 | 2.00 | 2.486 | 52% | FAIL ✓ |
| 8 | 0.125 | 1.00 | 0.50 | 3.243 | 68% | FAIL ✓ |
| 12 | 0.500 | 6.00 | 2.00 | 2.784 | 59% | FAIL ✓ |
| 12 | 0.250 | 3.00 | 1.00 | 4.820 | 102% | PASS ✓ |
| **16** | **0.500** | 8.00 | 2.00 | **2.773** | **58%** | **FAIL ✓** |
| 16 | 0.250 | 4.00 | 1.00 | 4.895 | 103% | PASS ✓ |
| 16 | 0.125 | 2.00 | 0.50 | 4.857 | 102% | PASS ✓ |

**20 of 20 correct.** Every predicted PASS measured 100-103%; every predicted FAIL
measured 41-86%. The corner named in advance is decisive:

    K =  8, Delta = 0.25  ->  100%        half the pilots, right spacing
    K = 16, Delta = 0.50  ->   58%        twice the pilots, wrong spacing

**Pilot count is not the variable. Pilot resolution is.** The design rule is

    (R1) Delta <= 2/P      comb finer than the sub-array beamwidth
    (R2) K*Delta >= 2      comb still spans the sweep

and with the baseline's hard-wired `Delta = 2/K` these collapse to **`K >= P`**,
i.e. **`K * N_TTD = Nt`** — the exchange law, now derived rather than fitted.

### 19.1 Correction to §14-§18: `K_min = P`, not `1.12 P`

K=8 passes at P=8. `frontier_ttd_pilot`'s ladder skipped K=7 and K=8 (it stepped
3, 6, 9, 12), so its `K_min = 9` was **ladder resolution, not a measurement**, and
`1.12 P` was the artifact of fitting four such ladder values. The law is

    K_min = P     <=>     K * N_TTD = Nt = 256

Confirmed at P=8. Still to confirm at P=4 (predict `K_min = 4`; the ladder
bracketed it in (3,5]) and P=16 (predict 16; bracketed in (12,18]).

### 19.2 BUG — the maxgap criterion had two counterexamples, and why

| K | `Delta` | gap/BW | rate |
|---|---|---|---|
| 4 | 0.125 | **0.6** | **41%** |
| 6 | 0.125 | **0.3** | **51%** |

Both read "pass" on a criterion I had reported as 25/25. The statistic was wrong,
not the criterion: it took `max(diff(sorted foci))` over the foci **alone**, so a
set that is dense but **narrow** — these two have comb span `K*Delta` = 0.5 and
0.75 against the 2 needed — piles its foci into a tight clump and reports a tiny
gap, because a difference taken only *between* foci cannot see the empty region
beyond the outermost one. Fixed in `cov_gap.m` by bounding the interval with the
region edges before differencing, so the distance from the extreme focus to the
edge counts as the hole it is. `decouple_theta`, `kspace_map`,
`diag_shared_mechanism` and `pilot_spacing_map` now all call it.

**The 25/25 claim must not be repeated until `pilot_spacing_map` is re-run with
the fix.** Earlier configurations all used the baseline comb (span 2) so their
numbers should be nearly unchanged, but that is an expectation, not a result.

### 19.3 What this opens — non-uniform pilot placement (`nonuniform_comb.m`)

A *uniform* comb provably cannot beat `K = P`: covering a span of 2 at resolution
`2/P` needs P points. But `map_focus_shift` showed the holes are **not uniformly
distributed** — they concentrate near boresight (P=4: `[-0.020, +0.137]`; P=8:
`[-0.140, +0.252]`). Where sharing opens no hole, fine resolution is wasted.

So a **non-uniform comb — dense where holes appear, sparse elsewhere — should
reach full coverage with `K < P`**, which no uniform comb can. And the selection
is free: each candidate offset's realized coverage comes from exactly the
`actual_focus` table the recalibrated lookup already computes offline, so this is
greedy set cover over precomputed sets, at zero run-time cost and zero extra
measurement.

- greedy `K < P` reaches ~100% while uniform at the same K does not → **pilot
  placement is a second algorithm-side lever**, "place them better" beats "use
  more", and `K >= P` is a property of uniform combs only;
- greedy never beats uniform → `K >= P` is fundamental to the beam family and
  `K * N_TTD = Nt` stands unconditionally, which is the cleaner paper.

---

## 20. CLOSED — one master curve, and non-uniform placement refuted

### 20.1 The gap statistic, re-run with the fix

Both counterexamples of §19.2 are gone: `K=4, Delta=0.125` and `K=6, Delta=0.125`
now read **128.7** beamwidths instead of 0.6 / 0.3, matching their 41% / 51% rates.
The criterion is **20/20 on this table with the corrected statistic**.

The 9 configurations of §16 and the 16 of §18 were scored with the buggy version.
All of them used the baseline comb (span 2), which has no narrow-clump case, so
their conclusions are unaffected — but **they have not been re-verified, and the
"25/25" figure must not be quoted** until they are.

### 20.2 The master curve — rate depends on `Delta*P/2` alone

Collecting every configuration that satisfies R2 (`K*Delta >= 2`), across K = 7
to 16 and three comb spacings, and plotting against the single ratio
**comb spacing / sub-array beamwidth**:

| `Delta*P/2` | configurations | rate |
|---|---|---|
| 0.50 | K=16, `Delta`=0.125 | 102%, 103% |
| 0.67 | K=12, `Delta`=0.167 | 101% |
| 1.00 | K=8/12/16, `Delta`=0.25 | 100%, 101%, 102%, 103% |
| 1.14 | K=7, `Delta`=0.286 | **98%** |
| 1.33 | K=6, `Delta`=0.333 | **73%** |
| 2.00 | K=4..16, `Delta`=0.50 | 52%, 55%, 56%, 58%, 59%, 59% |

**K itself has dropped out.** At `Delta*P/2 = 2` the rate is 52-59% for K = 4, 6,
8, 12 and 16 alike — a four-fold change in pilot count moves nothing. At
`Delta*P/2 = 1` it is 100-103% for K = 8, 12, 16. One variable explains the whole
table. **This is the paper's key figure**: a single curve of rate against
`Delta/(2/P)`, with every (K, `Delta`, P) point falling on it.

R1 is a **soft** threshold, as a beamwidth argument should be: 1.00 → 100%,
1.14 → 98%, 1.33 → 73%, 2.00 → ~55%. The knee is between 1.14 and 1.33.

### 20.3 `K_min` corrected once more: 7 at P=8

With the uniform baseline comb (`Delta = 2/K`, span 2) the transition sits between
K=6 (73%) and K=7 (98%), so at the 95% criterion **`K_min = 7 = 0.875 P`**, not 8
and not 9. The honest statement of the exchange law is

    K_min ~ P        <=>        K * N_TTD ~ Nt

with the constant measured at 0.875 for P=8 and the threshold soft. Earlier values
(1.12 P from §15, P from §19) were ladder artifacts of successively finer scans;
this one comes from a comb-spacing sweep and a rate curve rather than a bracket.

### 20.4 Non-uniform placement — REFUTED, and principled

`nonuniform_comb.m`, greedy set cover over 32 precomputed coverage sets:

| K | uniform | greedy |
|---|---|---|
| 5 | 64% | 66% |
| 6 | 71% | **76%** |
| 7 | **98%** | 77% |
| 8 | **100%** | 79% |
| 9 | **103%** | 87% |

Greedy wins by 2-5 points where both fail, and loses by 21-24 points everywhere it
matters. Its chosen offsets, sorted, are
`0.000 0.312 0.438 1.062 1.125 1.438 1.688 1.938` — **largest gap 0.624, two and a
half times the `2/P = 0.25` limit.** It broke exactly the property that §19 proved
decisive.

Two things went wrong, one fixable and one not:

- **fixable:** the coverage bins were 1-D in angle, while the requirement is 2-D
  over `(theta, alpha)`. Greedy optimized the wrong space.
- **not fixable:** even with correct bins, the required resolution is `2/P`, set by
  the sub-array beamwidth — a quantity that **does not vary with position**. A
  fixed budget of K points covering a span of 2 at a position-independent
  resolution is minimized by the **uniform** comb. Non-uniform placement can only
  help when the required resolution varies, and here it does not.

The holes concentrating near boresight (§17) does not change this: the holes live
in `sin(theta)`, the comb lives in offset space, and an offset shift translates the
whole strip pattern — so uniform spacing in offset is what produces uniform
resolution in coverage.

**Research question 2 is therefore answered in full, and negatively for the
"cheaper algorithm" route**: there is no placement trick. The pilots must be paid.

### 20.5 Status of the three research questions

| Q | verdict |
|---|---|
| 1. Does DD-BS's architecture/algorithm survive limited TTD? | **No.** 32 TTDs at baseline parameters → 42-43% of ideal. |
| 2. Design something that approaches ideal under limited TTD | **Yes, and the cost is now exactly characterized.** (i) recalibrated lookup — free, algorithm-only, and *necessary*: 85% of foci move by more than a grating step; (ii) pilot comb at `Delta <= 2/P`, i.e. `K ~ P`. No cheaper placement exists (§20.4). |
| 3. Abandon DD-BS for a new architecture | **Closed as posed.** The `theta_t` offset is load-bearing (§16.2) and the sweep cannot be re-centred away; the sharing limit is a sub-array-beamwidth property, not a parameter choice. |

---

## 21. The master curve holds across P; the criterion threshold is ~3, not 1; the OJ-COMS arm is unsafe

### 21.1 CONFIRMED — rate is a function of `P/K` alone, across a 4x range of P

`kmin_fine.m` was built so that two pairs of cells, measured at **different P**,
land on exactly the same `P/K`. The master curve says they must give the same rate.

| `P/K` | cell | rate | cell | rate | agreement |
|---|---|---|---|---|---|
| 1.333 | (P=4, K=3) | **76%** | (P=8, K=6) | **72%** | 4 pts |
| 1.143 | (P=8, K=7) | **99%** | (P=16, K=14) | **101%** | 2 pts |

Both agree inside Monte-Carlo noise (N_iter=30, ~±3 pts). **The master curve is
not an artifact of the P=8 sweep it was built from.** This is the single most
load-bearing validation in the project: it is what licenses quoting one curve
instead of one table per P.

### 21.2 `K_min = ceil(P/1.19)` — three for three

| P | K=… fails | K=… passes | `K_min` | `K_min/P` | `K*N_TTD` |
|---|---|---|---|---|---|
| 4 | 3 (76%) | **4** (100%) | 4 | 1.000 | 256 |
| 8 | 6 (72%) | **7** (99%) | 7 | 0.875 | 224 |
| 16 | 13 (78%) | **14** (101%) | 14 | 0.875 | 224 |

The crossing is bracketed between `P/K = 1.143` (passes) and `1.231` (fails), so
the threshold is **`P/K* = 1.19 +/- 0.04`**, tighter than §20's 1.14-1.33 bracket.
`K_min = ceil(P/1.19)` reproduces all three. My registered prediction
`ceil(P/1.14)` gave 4 / 8 / 15 and was **off by one at P=8 and P=16** — the
constant was slightly too small, which is exactly what a unit-step scan is for.

    K_min = ceil(P / 1.19)        K * N_TTD ~ 0.875 * Nt = 224   (P >= 8)

P=4's ratio of 1.000 is the ceiling rounding up from 3.36, not a different law.

### 21.3 Criterion threshold corrected: `gap <~ 3` beamwidths, not 1

The new cells put four points in the previously empty 1-5 band, and they move the
boundary:

| gap/BW | 0.3-0.9 | **1.4** | **2.9** | **3.2** | 5.2 | 8.9-10.1 | >20 |
|---|---|---|---|---|---|---|---|
| rate | 95-103% | **101%** | **99%** | **78%** | 80% | 72-89% | 41-64% |

Every pass has `gap <= 2.9`; every fail has `gap >= 3.2`. **Threshold ~3, still
zero misclassifications** — but the two nearest points are adjacent, so the
bracket is tight and the paper must quote `~3` with that bracket rather than a
round number. §16's "<= 1" was an artifact of having no data between 1 and 5.

Tally with the corrected statistic: `pilot_spacing_map` 20 + `kspace_map` 16 +
`kmin_fine` 8 = **44 evaluations, zero misclassifications**. Two of those
duplicate a cell measured by another script and agree to within one point
(K=6 s=1: 72% vs 73%; K=12 s=1: 101% vs 101%) — a free reproducibility check.

`kspace_map` and `diag_shared_mechanism` re-ran unchanged apart from the `span`
column (1.04 → 1.00, a cosmetic change in the `cov` formula), confirming the
§19.2 expectation that baseline-comb configurations were never affected by the bug.

### 21.4 Head-to-head — arm E dominates, at every equal-hardware, equal-overhead cell

`compare_ojcoms.m`, N_iter=30. Percentages are against arm A **at the same K**.

| N_TTD | K | `P/K` | A full | B shared/designed | C ojcoms/designed | D ojcoms/recal | **E shared/recal** |
|---|---|---|---|---|---|---|---|
| 64 | 3 | 1.33 | 4.748 | 2.031 (43%) | 1.588 (33%) | 1.344 (28%) | 3.527 (74%) |
| 64 | 4 | 1.00 | 4.831 | 2.829 (59%) | 1.887 (39%) | 1.428 (30%) | **4.800 (99%)** |
| 32 | 3 | 2.67 | 4.748 | 0.604 (13%) | 0.642 (14%) | 0.387 (8%) | 1.841 (39%) |
| 32 | 8 | 1.00 | 4.918 | 1.680 (34%) | 1.528 (31%) | 0.520 (11%) | **4.793 (97%)** |
| 16 | 8 | 2.00 | 4.918 | 0.917 (19%) | 0.565 (11%) | 0.292 (6%) | 3.124 (64%) |
| 16 | 16 | 1.00 | 4.968 | 1.462 (29%) | 0.934 (19%) | 0.330 (7%) | **4.838 (97%)** |

The headline that survives every caveat below:

> **16 TTDs (16x fewer than DD-BS) at K=16 deliver 4.838 bit/s/Hz — above the
> baseline's own 256-TTD, K=3 operating point of 4.748.**

`E - B` at fixed (N_TTD, K) isolates the **decision rule** on identical hardware:
29% → 97% at 16 TTDs. The recalibrated lookup carries most of the gain, and the
`P/K` law says when it is allowed to.

### 21.5 STOP — the `ojcoms` arm is unsafe and must not be published as it stands

Arm D is **worse than arm C at every cell** (7% vs 19%, 11% vs 31%, 30% vs 39%):
recalibration, which lifts every other configuration, *degrades* this one. That is
a signal about the model, not a result about their architecture.

Root cause in the code: with no fc compensation, the `ojcoms` intra-group residual
is `2*pi*f_m*(tau_g - tau_n)` — the **full** delay difference, not the
`(f_m - fc)` part. At the P=8 frontier that is ~42 carrier cycles of spread **at
every subcarrier including fc**, so the sub-array never coheres, the beam has no
focus to find, and `actual_focus` returns noise — which is worse to decode against
than a coherent-but-wrong designed table.

**A real AoSA cannot work that way, or OJ-COMS would not have published it.**
Almost certainly their phase shifters carry the `fc` part of the intra-group delay,
which makes their architecture **identical to the `shared` arm** — in which case
arms C and D are a strawman and must be deleted.

**Action required before any comparison goes in the paper**: check the OJ-COMS
text for whether the per-antenna PS term includes an `fc*tau` (or equivalent
sub-array-reference) compensation. Two outcomes:

- **it does** → arms C/D are deleted; their architecture *is* arm B/E's hardware,
  and the comparison becomes B vs E: same array, their decision rule vs ours,
  29% → 97%. The contribution is then **the algorithm and the design law**, not
  the array — a narrower but entirely defensible claim, and the one the data
  actually supports.
- **it does not** → arms C/D stand, but the paper must then explain how their
  published results are obtained despite the decoherence, which is a much heavier
  burden than the comparison is worth.

Until that check is done, **no number from arms C or D may be quoted anywhere.**
`ddbs_beam_arch.m`'s docstring now carries this warning at the `ojcoms` case.

---

## 22. PRIOR-ART FINDING #3 — the OJ-COMS paper already contains most of this project's results

Verified against the full text (IEEE OJ-COMS vol. 7, 2026, DOI 10.1109/OJCOMS.2026.3695965).

### 22.1 §21.5's hypothesis is REFUTED — the `ojcoms` arm is faithful

Their eqs (13)-(14), verbatim:

    (13) v^TTD_m(theta't, alpha't)|n_s = (1/sqrt(Nt)) exp(-j k_m [ n_s L d theta't - n_s^2 L^2 d^2 alpha't ])
    (14) v^PS_m (theta'p, alpha'p)|n   = (1/sqrt(Nt)) exp(-j k_c [ n   d   theta'p - n^2   d^2   alpha'p ])

The TTD carries the **full `k_m`** on the sub-array-centre lattice; the PS carries
**only** the DDBS phase parameters at `k_c`. **There is no `fc*tau` intra-group
compensation term.** `ddbs_beam_arch`'s `ojcoms` case reproduces their architecture
exactly, and the ~42 carrier cycles of intra-group spread are real, not a modelling
error. My §21.5 guess ("their PS almost certainly carries the fc part") was wrong.

They handle it not by compensating but by **accepting** it: the abstract states
they "characterize how delay sharing reshapes DDBS behavior by **restricting the
strong-gain region to a limited angular range/sector**", then cover the region with
sector shifting plus distance interleaving.

### 22.2 The overlap, with their equation numbers

| their result | equation / location | what I had called it |
|---|---|---|
| "three pilots ... are **not sufficient**" under AoSA | abstract | §12/§14 finding 1 |
| **LS effective focus** `theta_m = (kc/km) theta'p + rho_th theta't`, `alpha_m = (kc/km) alpha'p + rho_al alpha't`, `rho_th = L C1/S2`, `rho_al = L^2 C2/S4` | (29)-(31), Sec IV-B | **the recalibrated lookup — and in CLOSED FORM, no grid search** |
| Dirichlet intra-sub-array factor `|sin u / u|`, `u_3dB ~ 1.39`, `Delta_delta_strong = 5.56/L`, "scales as **O(1/L)**" | (43)-(44) | §19-§20 "comb finer than the sub-array beamwidth `2/P`" |
| **sector width** `W = 5.56 |S| / (L pi |theta't| xi_H^2)`, `S = theta'p + 2p/L` | (48), (53) | §16 equation (*) — the ratio law I declared refuted |
| **pilot count** `N_sec = ceil(theta_tot / W)`, hence `∝ L |theta't|/|theta'p|` | (54) | `K_min ∝ P` |
| `K_total = N_sec * K_alpha = 4*3 = 12` at `Nt=256, L=2` | Sec V, VII | my "K=12" headline |
| `N_TTD = Nt/L`; Table 2 against 2560 / 266 / 140 / 10 / 3 | Sec V-D, Table 2 | my hardware-vs-overhead table |
| **selection of L** as the hardware-overhead-performance knob | Sec V-E | my `(P, K)` frontier |

**That is essentially the whole contribution set assembled in §12-§21.** It must be
treated as prior art, cited, and not re-claimed in any form.

### 22.3 It also shows where I went wrong, twice, in the same place

- **§16's ratio equation was RIGHT.** `K >= Theta_req P |theta_t| / (4 eta |theta_p|)`
  is their (53)-(54) with `eta` supplied rigorously by the Dirichlet 3 dB solution
  `u_3dB = 1.39` rather than my hand-waved 0.443.
- **§18's "refutation" of it was an ARTIFACT of my own training design.**
  `decouple_theta` varied `theta_t` and `theta_p` while shifting pilots on a
  **uniform `theta_t` comb** — it never re-centred `theta'p` per pilot (their
  **sector shifting**, (62)-(64)) and never interleaved in `alpha` (their
  **distance interleaving**, (71)). With the sectors left stacked on top of each
  other, shrinking `theta_t` cannot help, so the ratio law looked dead. The same
  omission explains `nonuniform_comb`'s failure, which §20.4 correctly diagnosed as
  1-D angular binning against a 2-D `(theta, alpha)` requirement without connecting
  it to the missing second pilot dimension.
- Consistently, `kspace_map` scaled `theta_t` and `theta_p` **together**, holding
  their ratio fixed — so its finding that `s` does not matter is exactly what the
  ratio law predicts, not evidence against it.

### 22.4 What actually survives — one architectural difference, worth checking properly

Their PS (14) is programmed with `k_c(n d theta'p - n^2 d^2 alpha'p)`. The `shared`
arm programs the **same hardware** — one TTD per sub-array, one frequency-flat PS
per antenna — with that term **plus the exact per-antenna `2 pi fc tau_n`**. No
extra components, no extra resolution: only a different value written to a phase
shifter that already exists.

The intra-group residual changes from `2 pi f_m (tau_g - tau_n)` to
`2 pi (f_m - fc)(tau_g - tau_n)`. In their own derivation, (46) reads
`delta = pi xi_m theta't + 2 pi p / L`; under compensation the `pi xi_m` becomes
`pi (xi_m - 1) = pi (f_m - fc)/fc`, so at the band edge

    xi_H = 1.083   ->   xi_H - 1 = 0.083          a factor of ~13 at B/fc = 1/6

Through (52)-(54) that is a **~13x wider sector and ~13x fewer sector shifts**.
Their `L=2` operating point needs `N_sec=4, K_alpha=3, K_total=12`; the same
argument would put the compensated architecture at `N_sec=1`, i.e. `K_total = 3`.
At `L=8` their formula gives roughly `N_sec ~ 16, K_total ~ 48`, against a measured
`K_min = 7` for the compensated arm.

Measured, at equal `(N_TTD, K)` (§21.4): 39% → 99%, 31% → 97%, 19% → 97%.

**But those numbers do NOT compare against their scheme.** Arm C is their
architecture driven by *my* uniform-comb training, not their sector-shifted,
distance-interleaved, LS-refined procedure. Quoting 19% vs 97% against their paper
would be a false comparison and must not be done.

### 22.5 Required next step, and nothing may be claimed before it

Implement **their full scheme** as the baseline — (29)-(31) LS focusing, (53)-(54)
sector shifting, (62)-(64) sector parameters, (71) distance interleaving — at their
operating point (`Nt=256, L=2, N_sec=4, K_alpha=3, K_total=12`), and reproduce their
reported performance. Only then add the fc-compensated PS and compare at equal
`(N_TTD, K_total)`. Two outcomes:

- **compensation reduces `K_total` substantially at equal `N_TTD`** (the ~13x
  argument, testable directly through their (53)) → a genuine, narrow, defensible
  contribution: *the same AoSA hardware, reprogrammed, needs far fewer pilots* —
  with the delay-phase split itself cited to Dai/Tan TWC 2022 (§13) and the
  sector/coverage framework cited to this paper;
- **it does not** → this project has no publishable increment over OJ-COMS 2026,
  and that conclusion should be reached quickly rather than slowly.

Everything in §12-§21 stays valid as *measurement*; what changes is that almost
none of it is *novel*, and the write-up must be rebuilt around 22.4 alone.

---

## 23. Gate FAILED, diagnosed, and Table 3 obtained

`ojcoms_baseline.m` stage 1, first attempt:

| L | K | O-cf | O-af | their Fig. 7 |
|---|---|---|---|---|
| 1 | 3 | 0.212 | 0.635 | ~6.5 |
| 2 | 12 | 0.258 | 1.503 | ~5.0 |

Off by an order of magnitude — the reconstruction was wrong, exactly the case the
gate exists to catch. Two parts of it were right and two were bugs.

### 23.1 What was right

- **`rho_theta = 1.0000, rho_alpha = 0.9999` at L=2 is correct**, not a degenerate
  case. Analytically `C1 = 4*sum_k (k-63.5)^2` for L=2, giving
  `rho_theta = 2*699008/1398080 = 0.99995`. The LS coefficients only depart from 1
  at larger L.
- **`alpha't = -0.5338` against their published `-0.533`.** (69)-(70) and `U` are
  implemented correctly.

### 23.2 Bug 1 — `S` solved from the wrong branch

I set `S` by boundary focusing (`theta_1 = -1, theta_M = +1` on one strip),
getting `S = -11.92`. Their text says `S_m^(s)` is chosen **near the upper bound**
of (59)-(60), i.e. `S ~ +S_bound ~ 34.85`. Feeding `S = 35.53` with `theta_M = 1`
into their (63) gives

    theta't = (1 - (fc/fH)*35.53)/rho_theta = 1 - 0.92308*35.53 = -31.797

— **their published value exactly.** That is the recipe. With `gamma = 0.9`,
`S = S_bound = 34.85` gives `-31.17`, 2% off, which is a reproduction rather than
a fit. The physical reading: they run a **fast** sweep with many alias wraps
inside the band (`S ~ 35` moves the raw focus by ~5.9 in `theta` across the band,
i.e. ~6 strips per pilot), not one slow strip across the range.

### 23.3 Bug 2 — the closed-form lookup dropped the alias term

(29)-(30) give the **principal branch only**. With `theta't ~ -32` and
`theta'p ~ 28` the principal focus sits near `-3.4`, far outside `[-1,1]`, so the
O-cf arm was pointing nowhere — hence 0.212. The realized focus needs the alias
term of (39)-(40), and **their own (43)/(46) say which alias wins**: the
intra-sub-array factor is `|sin(u)/u|` with `u = L*delta/2` and
`delta = pi*xi_m*theta't + 2*pi*p/L`, so the strongest branch is the `p` that
maximises it. That is now implemented, restricted to candidates inside `[-1,1]`.

### 23.4 Table 3, obtained — the baseline is now THEIR configuration

| sector | `theta't` | `theta'p` | `gamma` | `alpha't` | `alpha'p` |
|---|---|---|---|---|---|
| 1 | -31.797 | 28.40 | 0.9 | -0.533 | 0.158 |
| 2 | -33.840 | 31.93 | 0.9 | -0.533 | 0.158 |
| 3 | -33.545 | 31.61 | 0.9 | -0.533 | 0.158 |
| 4 | -32.102 | 28.51 | 0.9 | -0.533 | 0.158 |

Used verbatim for L=2, so the comparison baseline is their published
configuration and not a reconstruction of it.

### 23.5 Two things the paper does not determine, flagged rather than guessed

1. **The split of `S` into `(theta'p, p_M)`.** Their "convenient implementation"
   `p_M = floor((L/2) S_bound) = 34` gives `theta'p = 0.85`, not the 28.4 of
   Table 3, which corresponds to `p_M = 7`. Only the composite `S` is recoverable
   from the text. Table 3 settles it for L=2; for other L the analytic branch
   carries their L=2 split forward, which is an assumption.
2. **What separates the `K_alpha = 3` distance pilots inside one sector.** Table 3
   lists **one** `(alpha't, alpha'p)` per sector and the text says they are "fixed
   across all sectors" — taken literally the three pilots in a sector are
   identical. Interpreted here as interleaving `alpha't` across the (70) interval.
   `Kalpha = 1` reproduces Table 3 literally and is now also run, so the gate
   reports both 4-pilot and 12-pilot results.

Until the re-run passes, **no comparison number may be quoted**, and §21.4's
39%/31%/19% column stays unusable.

---

## 24. Gate: L=1 passes; L=2 exposes a contradiction inside their own paper

Re-run after the §23 fixes, with Table 3 used verbatim:

| L | K | O-cf | O-af | their Fig. 7 |
|---|---|---|---|---|
| 1 | 3 | **5.857** | **5.938** | ~6.5 |
| 2 | 4 | 0.831 | 0.337 | ~5.0 |
| 2 | 12 | 1.109 | 0.765 | ~5.0 |

The unit test now reproduces `theta't = -31.797, theta'p = 28.400` exactly, and
`alpha't = -0.5338` against their `-0.533`.

### 24.1 The ungrouped path is validated

L=1, K=3 gives **5.86-5.94 against their ~6.5, i.e. 90%**. Channel model, rate
metric, operating point (`r` in [5,200] m, SNR 20 dB) and DDBS beam generation
are all consistent with their setup. The residual 10% is plausibly the user
distribution, which their text does not specify (uniform in `r` assumed here).

**So the failure is specific to the grouped path**, which is what matters.

### 24.2 Table 3's `alpha'p` contradicts their own equation (69)

With Table 3's `alpha'p = 0.158` and `alpha't = -0.533`, the LS focus of (66)
sweeps

    alpha_m in [-0.3872, -0.3606]

— the **entire pilot set points outside the target interval `[0.0025, 0.1]`**, so
no user is reachable and the measured rate says nothing about the architecture.

Their own (69) gives

    U = (alpha_max - alpha_min)/(fc/fL - fc/fH) = 0.5809   ->   alpha_m in [+0.0033, +0.1008]

which is exactly the target interval — by construction, since that is what (69)
is for. Three independent confirmations that **0.158 is a typo**:

1. (69), their own equation, gives 0.5809;
2. 0.5809 sweeps precisely `[alpha_min, alpha_max]` while 0.158 sweeps nothing usable;
3. 0.5809 is also the DD-BS baseline's own `alpha_p = 1859/3200 = 0.5809`, and
   their `alpha't = -0.533` likewise matches the baseline's `-427/800 = -0.534`.

So the pair is `(alpha't, alpha'p) = (-0.533, 0.581)` and Table 3 mis-set the
second. Using (69) in preference to a table entry that the paper's own equation
contradicts is a correction, not a tuning step — and both are now run, so the
out-of-range behaviour of the table value is documented rather than hidden.

`ojcoms_baseline.m` now runs the gate twice and prints the realized `alpha` sweep
next to each rate. **If the eq-(69) block brings L=2 near 5.0, the baseline is
validated and stage 2 may proceed.** If it does not, the remaining suspect is the
`(theta'p, p_M)` split of §23.5, and the honest conclusion would be that their
scheme is not reproducible from the published text.
