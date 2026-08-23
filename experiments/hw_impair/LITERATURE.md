# Literature map — near-field wideband beam training with limited TTD hardware

**Caveat on the tiering.** The 中科院分区 below are my best estimate from the
journal, not verified against the current 分区表 — it changes yearly and by
大类 (工程技术 vs 计算机科学). **Check every one against the latest table before
citing a tier in a proposal.** Where I have not confirmed the venue of an arXiv
preprint, it says so.

---

## A. The direct lineage — must all be cited

| # | work | venue | tier (est.) | why it matters here |
|---|---|---|---|---|
| A1 | **Near-field wideband beam training based on distance-dependent beam split** (arXiv 2406.07989) | IEEE **TWC** | **1区** | the baseline. K=3, one TTD per antenna |
| A2 | Cui & Dai, **Near-field rainbow: wideband beam training for XL-MIMO** (10.1109/TWC.2022.3222198) | IEEE **TWC** 2023 | **1区** | 10 pilots, per-antenna TTD. OJ-COMS's `[22]` |
| A3 | Dai, Tan, Chen, Poor, **Delay-phase precoding for wideband THz massive MIMO** (arXiv 2102.05211) | IEEE **TWC** 21(9) 2022 | **1区** | **the TD/PS split itself — prior art for contribution A (THEORY §13)**. Also `[31]` in A1's own reference list |
| A4 | Nguyen & Kim, **Joint delay-phase precoding under true-time-delay constraints in wideband sub-THz hybrid massive MIMO** (10.1109/TCOMM.2024.3402616; arXiv 2212.07484 is the preprint) | IEEE **TCOM** 72(10), Oct 2024 | **2区** (工程技术; check the year) | **READ — see THEORY §25.** Same AoSA structure; joint TTD **and** PS optimization under bounded per-TTD delay, explicitly against fixing the PS; closed-form **minimum-TTD count** for a target array gain, growing **linearly with B**. Far-field, squint *compensation*, no training — so it does not block the claim, but it owns "use the PS properly" as a principle |
| A5 | Najjar, El-Absi, Kaiser (TTD resolution constraint) | IEEE **TWC** 23(7) 2024 | **1区** | derives `phi = 2 pi fc tau` for the finite-resolution motivation — prior art (THEORY §13) |
| A6 | Qaid, Nasir, Al-Ahmadi, Liu, **A hardware-efficient hybrid beamforming architecture…** (10.1109/OJCOMS.2026.3695965) | IEEE **OJ-COMS** v7 2026 | **3区** | **the paper this project must beat.** AoSA grouped TTD, sector shifting, LS focusing, MF refinement — see THEORY §22 |

## B. Same problem, adjacent architecture — the competitive field

| # | work | venue | tier (est.) | note |
|---|---|---|---|---|
| B1 | Wadaskar et al., **Fast 3D beam training with TTD arrays in wideband mmWave** | IEEE **TWC** 24, pp.5146-5162, Jun 2025 | **1区** | far-field 3D, same TTD-rainbow family |
| B2 | Boljanovic, Yan, … Cabric, **Fast beam training with TTD arrays** (IEEE Xplore 9349090) | IEEE **TCAS-I** 2021 | 2区 | the origin of far-field rainbow training |
| B3 | **Near-field wideband beamforming for extremely large antenna arrays** (10.1109/TWC.2024.3398770) | IEEE **TWC** 2024 | **1区** | near-field wideband precoding, not training |
| B4 | **Beamfocusing optimization for near-field wideband multi-user** (arXiv 2306.16861) | venue unconfirmed | — | OJ-COMS cites it for **antenna-group** beamfocusing |
| B5 | **Energy-efficient dynamic-subarray with fixed TTD for THz wideband hybrid beamforming** (arXiv 2202.02965) | venue unconfirmed | — | **fixed** (non-tunable) TTD + dynamic sub-arrays — the cost axis we also touch |
| B6 | **Efficient hybrid near- and far-field beam training for XL-MIMO** (Xplore 10643414) | IEEE 2024 | — | hybrid-field training, no TTD |
| B7 | **Wideband near-field channel covariance estimation … in the face of beam split** (Xplore 10709906) | IEEE 2024 | — | estimation rather than training |
| B8 | **Time-delay-aided fast beam training for near-field** (10.26599/TST.2025.9010112) | Tsinghua Sci & Tech 2025 | 2-3区 | TTD-based 2-phase frequency scanning |

## C. JPTA — the industrial line on frequency-dependent analog beamforming

| # | work | venue | tier (est.) | note |
|---|---|---|---|---|
| C1 | **Joint phase-time arrays: a paradigm for frequency-dependent analog beamforming in 6G** (arXiv 2312.11682) | IEEE **Access** | 3区 | Samsung. The name OJ-COMS models its PS on (`[23]`) |
| C2 | **Hybrid near/far-field frequency-dependent beamforming via JPTA** (arXiv 2501.15207) | unconfirmed, 2025 | — | near-field JPTA, single RF chain |
| C3 | **JPTA: opportunities, challenges and system design considerations** (arXiv 2412.01714) | unconfirmed | — | design-constraint survey |
| C4 | **Beamforming with joint phase and time array: system design, prototyping** (arXiv 2502.00139, Xplore 10942922) | IEEE conf 2025 | — | **hardware prototype** — useful for the realizable-delay-range argument of §13 |

## D. Surveys, for the intro

- **Near-field communications: a tutorial review** (arXiv 2305.17751)
- **Cross far- and near-field beam management in mmWave and THz MIMO** (arXiv 2504.18855)
- IEEE ComSoc **Best Readings in Near-Field MIMO Communications**

---

## What this map says about the remaining contribution

The AoSA arrangement we call the `shared` arm — **TTD on the sub-array common
leg, per-antenna phase shifters absorbing the intra-sub-array delay at `fc`** —
is **standard phased-array practice**, described in vendor application notes and
in the wideband-AoSA literature (and it is what A3/A5 formalize). **It must not
be claimed as a new architecture.**

What is not in the map:

1. **Nobody applies it to DDBS *training* under sub-array sharing.** A6 is the only
   grouped-TTD DDBS training paper, and its (14) deliberately does **not**
   compensate at `fc`.
2. **A6 reports the resulting ceiling as a property of grouping**: their Fig. 7
   gives ~6.5 / 5.0 / 3.7 / 1.9 / 1.3 / 0.5 bit/s/Hz at `L` = 1/2/3/5/7/9 even
   with the pilot budget their own (73) requires, and they conclude `L=2` is the
   favourable operating point. If the compensated arrangement removes that
   ceiling, then the ceiling was **an artifact of their phase-shifter design, not
   of delay sharing** — and that is the claim worth making.
3. The pilot law `K_min ~ P/1.19` and the offline `maxgap` criterion (§20-§21)
   are ours as *measurements*, but their §V already gives an analytical pilot
   count (their (53)-(54)); ours would have to be positioned as **the law for the
   compensated architecture**, which theirs does not cover.

**Priority read: A4** (arXiv 2212.07484). It is the closest published work to the
surviving claim, and if it already treats sub-array-shared TTD with `fc`
compensation under a delay-range constraint, the contribution narrows again.


---

## Update after reading A4 (THEORY §25)

A4 is **IEEE TCOM 72(10), Oct 2024** — a stronger venue than A6, and it owns the
joint-TTD-PS principle plus a min-TTD-count law. The surviving contribution is
therefore an **application** of a known principle to a problem where it had not
been applied, not a discovery:

- **inherited**: the AoSA architecture (vendor practice); the delay-phase split
  (A3, TWC 2022); the joint design principle and min-TTD law (A4, TCOM 2024); the
  sectoring/coverage framework and the grouped-DDBS training problem (A6, 2026).
- **remaining**: one architectural substitution inside A6's scheme, plus the
  measured pilot law `K_min ~ P/1.19` and the offline `maxgap` design test —
  with the **bandwidth-independence** of our law standing in deliberate contrast
  to A4's `B`-linear TTD law.

**Realistic target: 2-3区** (OJ-COMS-tier, or a solid 2区 such as TVT/TCOM if the
execution is strong). Scope it that way from the start.

---

## E. Beam-squint-assisted LOCALIZATION / ISAC — read, and they constrain us

| # | work | venue | tier (est.) | bearing on this project |
|---|---|---|---|---|
| E1 | Luo & Gao, **Beam squint assisted user localization in near-field ISAC systems** | IEEE **TWC** 23(5), May 2024 | **1区** | **Must cite.** They *"derive the trajectory equation for near-field beam squint points and design a way to control such trajectory"*, then use the frequency-domain squint to localize with reduced sweeping overhead — the same physics as DDBS strips, applied to sensing. Also `[44]` in Nguyen & Kim TCOM 2024. **Refinement idea:** multi-carrier **phase difference** between max-power subcarriers across several sweeps → distance RMSE ~0.10 m at 15 dB. |
| E2 | **Beam-squint assisted joint angle-distance (JAD) localization for near-field communications** (DOI 10.1109/TVT.2026.3706538, accepted) | IEEE **TVT** 2026 | 2区 | **Must cite.** Coarse-to-fine two-stage: stage 1 = coarse joint `(theta, r)` from the **power-spectrum peak** (this is exactly our decision rule); stage 2 = **near-field improved MUSIC** searching locally around it. Explicitly motivated by the **error propagation** of two-step angle-then-distance estimators. |

### What they mean for the claim — the estimator is not available to us

Both papers already occupy the **estimator** territory this project might have
moved into:

- E1 owns "control the near-field squint trajectory and read location out of the
  frequency domain";
- E2 owns "coarse peak, then joint high-resolution refinement, avoiding
  sequential error propagation".

So **novelty cannot live in the decision rule.** Adding MUSIC or a phase-based
refinement would be re-doing E2/E1, on top of OJ-COMS's own MF refinement (their
Sec VI). The correct use is the opposite: **borrow the best available refinement,
apply it to BOTH arms, and cite them.** That turns it from a novelty risk into a
robustness result — a reviewer cannot then say the architecture gap would vanish
under a stronger estimator, because it was measured with one.

### The error-propagation question, answered

The intuition was right but it does not apply to this baseline:

- **DD-BS and OJ-COMS decide with a single `argmax` over subcarriers**, then read a
  joint `(theta, alpha)` from one table. There is no sequential stage to
  propagate into, so "we fix the baseline's error propagation" would be **factually
  wrong** and a reviewer would see it from the algorithm listing.
- Error propagation is real in **two-step angle-then-distance** estimators — and
  E2 already published the fix. It cannot be claimed here.
- What *can* be said, in one sentence citing E2: our decision rule is joint by
  construction, so it is structurally free of the effect E2 addresses. That is a
  remark, not a contribution.

### Two things worth borrowing, with their costs

1. **Phase-difference distance refinement (E1).** Our pipeline uses `|sum(y)|^2`
   and **throws the phase away**. Using the phase across max-power subcarriers
   could sharpen `alpha` with **no extra pilots** — and `alpha` is exactly where
   L=16 is weakest (§27, 89% at `L/N_sec = 1.33`). Cheapest transfer, highest fit.
2. **Coarse-to-fine local refinement (E2).** Maps 1:1 onto our pipeline: argmax +
   recalibrated table = their stage 1. **But their stage 2 needs a covariance from
   spatial snapshots**, and a single-antenna UE in beam training receives one
   scalar per (subcarrier, pilot) — not an array snapshot. Whether an equivalent
   covariance can be formed from the `M x K` measurements must be checked before
   adopting it; do not assume it transfers.
