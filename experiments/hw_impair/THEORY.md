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
