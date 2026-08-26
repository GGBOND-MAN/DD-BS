# `theory/` — derivation code for THEORY.md §40

A NumPy port of `ddbs_beam_arch.m` and `actual_focus.m`, used to derive and test
the coverage law. **No channel model and no baseline MATLAB code are involved** —
both ported functions are self-contained, so the focus geometry is reproducible
anywhere. Run from inside this directory (the scripts `exec` each other by
relative path):

| script | what it establishes | §40 |
|---|---|---|
| `verify_closed_form.py` | the ported `actual_focus` matches the closed-form locus to **1.3e-4** at `L=1`, once the alias period `2*fc/f_m` is included | 40.1 |
| `alias_orders.py` | which grating orders `j` are actually visible: `{-2,-1,0}` at `s=1`, `{11}` at `s=1/4`, `{13}` at `s=1/8`; and the `L=8` focus shift | 40.2-40.3 |
| `cover.py` | the coverage statistic `P(g>0.9)` on the twelve §38 cells, against the measured MATLAB rates (Spearman 0.907) | 40.4 |
| `collapse.py` | the pass/fail test of the derived law over `N_sec` x `s`, **12/12 with zero fitted parameters** | 40.5 |
| `confirm.py` | the forward prediction: at `s=1/8` sectors alone need `N_sec ~ 42` | 40.5 |
| `diag_ka.py` | the `K_alpha` gain is near-uniform across the user `alpha` deciles — not a density effect | 40.6 |
| `diag_clamp.py` | **40%** of the interleaved pilot's subcarriers are pinned to `alpha_min` over 2/3 of its angular slice | 40.6 |
| `diag_recal.py` | that mechanism needs the recalibrated lookup: `0.935` with it, `0.456` without | 40.6 |

`ddbs.py` is the port itself (`beam`, `focus_table`, `alg1`, `rho`).

The coverage statistic is
`g(user) = max over (pilot, subcarrier) of the narrowband focusing gain at fc`
between the user's `(theta, alpha)` and the **recalibrated** focus-table entry —
i.e. the geometric quality of the best beam the training stage can hand the
serving stage. It is a proxy, not a rate: it validates against the twelve
measured rates at Spearman **0.907**, with one inversion in the marginal pair
(`P=0.925 -> 98%` vs `P=0.936 -> 92%`). The direct end-to-end confirmation is
`hw_impair/kmin_theory.m`, which must be run in MATLAB against the real channel.
