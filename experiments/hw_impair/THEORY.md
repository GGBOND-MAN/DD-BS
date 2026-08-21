# TD resolution requirement for DDBS beam training — derivation

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
  ideal values: count x precision is a 2-D hardware-cost surface, and the two
  contributions compose.

---

## 7. Second axis: TD **count** — the Kronecker two-stage TTD

### Why naive sub-array sharing is expensive for DDBS
Sharing one TTD across `P` antennas leaves a residual `2*pi*(f_m-fc)*(tau_g - tau_n)`.
The intra-group delay spread is `P*d*|theta_t|/c`, and `|theta_t| ~ 33` for DDBS,
so even `P = 2` gives ~1.1 ns of spread — several radians at `B/2`. Measured
focusing gain: `P = 2 -> 0.668`, `P = 4 -> 0.402`, `P = 8 -> 0.238`.

### The structure that fixes it
Split the antenna index exactly as `n = n_g + n_q` (`n_g` the sub-array offset,
`n_q` the position within it). Then

    tau_n = (n*d*theta_t - n^2*d^2*alpha_t)/c
          = [A(n_g) + C(n_q)]  -  2*n_g*n_q*d^2*alpha_t/c ,

i.e. the **linear term is exactly separable** and only the quadratic **cross term**
resists. Since the DDBS delay is linear-dominated (140 ns vs 0.72 ns here — the
same large `theta_t` that makes the hardware hard is what makes this factorization
near-exact), a two-stage bank of

    Nt/P coarse TTDs  +  P fine TTDs (reused by every sub-array)  =  Nt/P + P

reproduces the delay profile. The count is minimised at `P = sqrt(Nt)`, giving
**`2*sqrt(Nt)` TTDs** (32 instead of 256 at `Nt = 256`, an 8x reduction).

### Design rule for P
The residual is bounded by `d_tau <= Nt*P*d^2*|alpha_t| / (2*c)`, so the worst-case
phase error is

    **eps_max = pi * B * Nt * P * d^2 * |alpha_t| / (2*c)**

— it grows **linearly in P**, while the TTD count falls as `Nt/P + P`. Choose
`P = min( sqrt(Nt), 2*c*eps_budget/(pi*B*Nt*d^2*|alpha_t|) )`. At the paper's
settings both terms land near `P = 16` (`eps_max = 1.4 rad`, gain 0.941).

### Measured surface (ideal TD, hybrid PS; ideal full-TTD = 0.975)

| #TTD | 128 | 64 | 32 |
|---|---|---|---|
| generic sharing | 0.668 | 0.402 | 0.238 |
| **Kronecker** | **0.975** | **0.973** | **0.941** |

~4x higher gain at equal TD count. And the two contributions **compose**:
Kronecker at 32 TTDs still gives **0.938 at 12-bit** TD — ~8x fewer TTDs *and*
~4 fewer bits each, versus 256 full-resolution (~14-bit) TTDs.
