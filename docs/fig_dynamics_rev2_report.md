# fig_dynamics rev2 — report

Branch `step5-kd-fits`. Data: corrected master export, stamp **`20260824_221720`**.

**PR #22 was UNMERGED at the start of this session, contrary to the handoff.** The rev2
precondition (`step5-kd-fits` contains `fix-apg-export`) therefore failed and was reported
before any edit. #22 was then merged (`9c883a8`), `step5-kd-fits` rebased onto main (clean,
7 commits replayed), and rev2 proceeded from that base.

---

## 1. What changed in each figure

**`fig_tstop`** — rebuilt to the KD Fig. 1b construction. 1×3 → **one panel**, all three
geometries overlaid, x from 0 to 300 cm/s. All raw trial markers removed; one point per
(geometry, height) cell, mean v₀ and mean t_stop_s with ±1 SD bars in both axes. Added the
dotted black ODE curve from the locked (k₀, c₀), and the dashed T_c / V_c reference lines
with numeric annotations. `tstop5` variants appear nowhere.

**`fig_scaling_models`** — unchanged except legend text: each entry now carries its equation
and R². The KD entry uses the **citation form**, not the equation: at single-column width the
full `m z̈ = mg − kz − c ż²` entry overflowed the axes, and the brief's fallback citation form
overflowed further (LaTeX `~` and `---` also rendered literally). Final text is the compact
`KD Eq. (1), no free parameters, R² = 0.639`; the full equation belongs in the caption.

**`fig_force_decomposition`** — full rebuild on KD Fig. 3 logic. Top row = three panels
(a+g vs v) with four depth colours, dotted fixed-slope parabolas, open intercept symbols at
v = 0, and a **zoomed** v² inset. Bottom = one combined panel (a+g − v²/d₁ vs z) with cell
curves, the 12 intercepts with CI and depth-window bars (±0.015 cm geometry offset), and the
two candidate F(z)/m forms. `c_bins` and `parallel p` removed from the panels to the report.
Dashed vertical lines removed.

**`fig_parameter_ladder`** — panel (c) and everything feeding it deleted; ΔBIC and selection
fractions are text-only (below). Panels (a)/(b) keep one summary marker per geometry at the
median of its 18 cell fits, with 2.5–97.5 percentile error bars; scatter clouds removed. Units
corrected to `k (10⁵ g/s²)` and `c = m/d₁ (g/cm)`. `SHOW_M0_REF = false`.

## 2. Reuse, and what was added

Reused: the fixed-depth bin extraction rule of `step5_depth_term_study`; the RK4 KD integrator
and root-finding closure of `step5_closure_diagnostics` (reproduced verbatim with the source
named at each site — MATLAB local functions are not importable across files); the persisted
cluster resamples of `step5_geometry_ladder`; the ensemble-median + support-gating construction
of `step5_fleet_fits`; the geometry identity block of `fig_scaling.m` (lines 101–106, verified
identical); `net_accel` via the exporter.

Added display-support computations only: per-depth intercepts at a fixed slope, the pooled
12-group common slope, ODE integration for the t_stop curve, resample percentiles, and the
KD-prediction R². **No fit parameter changed.** (k₀, c₀) remain the manuscript's law parameters.

Two brief items could not be honoured as written and were substituted, both stated on screen:

- `fig_scaling.m` defines **no** per-cell aggregation convention (it plots raw markers only),
  so `fig_tstop` uses the brief's first-named option — per-height mean with ±1 SD.
- The 50 % support rule attributed to `fig_kinematics` is actually `step5_fleet_fits`'
  construction; `fig_kinematics` gates on an absolute `MinReplicates = 3`. The 50 % rule was
  implemented as specified.

## 3. Reruns on the corrected export

Re-running the exporter from the rebased branch produced stamp `20260825_002752`, whose
trials CSV, series CSV and dictionary are all **byte-identical** to `20260824_221720`.
**`20260824_221720` remains authoritative**; the redundant duplicate was left in place.

**c_bins and parallel p — corrected vs pre-correction.** Reproduced exactly, as a constant
a+g offset cannot change a slope:

| geometry | c_bins (corrected) | c_bins (pre) | p (corrected) | p (pre) |
|---|---|---|---|---|
| Tight | 14.14 | 14.1 | 0.0725 | 0.073 |
| Default | 12.32 | 12.3 | 0.3393 | 0.339 |
| Wide | 13.94 | 13.9 | 0.0001 | 0.000 |

**c_a (new, primary for Fig. 4).** Pooled common slope over all 12 groups:
**c_a = 13.392 g/cm [11.901, 14.263]**, d₁ = 4.8536 cm.
Pooled parallelism (do the 12 slopes differ?): **F = 3.436, p = 9.36 × 10⁻⁵** — they do.
c_a is **−15.8 %** vs c₀ = 15.898, and **the c_a CI does not contain c₀**.

**Intercept table F(z_i)/m [cluster CI], `D1_SOURCE = 'panel'` (d₁ = 4.8536 cm):**

| geometry | z = 0.3 | z = 0.5 | z = 0.7 | z = 0.9 |
|---|---|---|---|---|
| Tight | 1705 [1274, 2137] | 1835 [1387, 2244] | 1775 [1353, 2201] | 1650 [1101, 2295] |
| Default | 2136 [1575, 2729] | 2866 [2346, 3367] | 2814 [2367, 3234] | 2854 [2534, 3200] |
| Wide | 2600 [2092, 3066] | 3124 [2512, 3693] | 3000 [2597, 3426] | 2697 [2261, 3133] |

n = 170–177 trials and 18 cells in every bin. The median-based intercept differs from OLS by
> 10 % in 6 of 12 bins (all of Default, and Wide at 0.3/0.7/0.9) — the bin distributions are
skewed, so the mean intercept sits above the median. Trial-level outliers beyond 3× IQR,
retained in all fits: Tight 0, Default 0, Wide 4.

**Collapse-quality Spearman(y(z), cell median v₀):**

| z | d₁ = m/c_a | d₁ = m/c₀ |
|---|---|---|
| 0.3 | +0.216 | −0.167 |
| 0.5 | +0.496 | +0.098 |
| 0.7 | +0.422 | −0.164 |
| 0.9 | −0.065 | −0.457 |

Neither d₁ removes the v₀ ordering: c_a leaves a positive residual ordering at 0.3–0.7 and c₀
leaves a negative one at 0.3/0.7/0.9. Descriptive only; the primary d₁ was not changed on this
basis.

**4B level shift.** The rev1 row-2 panel already used the corrected column, so the export fix
contributes ≈ 0. The shift is entirely the d₁ change, (c₀ − c_a)v²/m ≈ **+940 cm/s² at
150 cm/s** and **+2610 at 250 cm/s** — arithmetic, not a new departure.

## 4. Values

- d₀ = 0.0441, a = 0.07106, n = 0.6059.
- R²: fixed-2/3 **0.7051**, free-n **0.7052**, KD **prediction score 0.6385** (not a fit).
- T_c / V_c by choice of L_c: **D_eq = 1.644 cm → 0.0410 s / 40.1 cm/s** (default);
  √A_bare = 1.457 → 0.0386 / 37.8; ℓ₀ = 0.862 → 0.0297 / 29.1; d_s = 0.751 → 0.0277 / 27.1;
  d₁ = 4.089 → 0.0646 / 63.3.
- ODE t_stop at v₀ = 5 cm/s = **0.0529 s**; drag-free harmonic limit π√(m/k₀) = 0.0579 s.
- KD d(v₀→0) = **0.6323 cm** with c₀; **0.6648 cm** as c → 0; drag-free 2mg/k₀ = 0.6648 cm.
- t_stop fallback fraction: Tight 98.9 %, Default 98.3 %, Wide 98.9 %.
- Cell median vs pooled M0: k **+6.3 % … +27.1 %**; c **−18.8 % … −10.0 %**.
- Ladder (text-only): wRMSE {5000.9, 4995.4, 4998.2, 4978.7}; BIC {8938.7, 8950.1, 8950.7,
  8959.1}; ΔBIC vs M0 {0, 11.4, 12.0, 20.4}; selection {97.4, 0.2, 2.0, 0.4} %.
- 4B max plotted depth: Tight 3.05, Default 3.08, Wide 3.15 cm.

## 5. Caveats affecting manuscript use

1. **c_a vs c₀ — a real cross-method gap.** c_a = 13.392 [11.901, 14.263] vs c₀ = 15.898:
   **the CI excludes c₀** (−15.8 %). The likely explanation is that the windowed derivative
   attenuates a+g most at high v, biasing the acceleration-space slope low — the collapse
   Spearman table is consistent with over-correction under c_a — but this is not asserted.
2. **Wide non-parallelism (p = 0.0001)**, and the pooled 12-group test also rejects
   (p = 9.4 × 10⁻⁵). A single common slope is an approximation; the departure is absorbed into
   the intercepts, which is why they must not be read as k.
3. **Truncation.** Cell curves run from the first to last grid depth with ≥ 50 % of the cell's
   trials contributing, and never past the cell's median d_final. The rev1 abrupt ends near
   2.7–3.0 cm were per-cell depth range plus support gating; the rule is now uniform.
4. **L_c choice.** D_eq = 2√(A_bare/π) is the literal D_b analogue and is identical for all
   three feet; the four alternatives are tabulated above and L_c is a single named parameter.
5. **ℓ₀ mapping.** ℓ₀ = 0.862 reproduces as `b^{3/2}·√(2g)` from the pure-2/3 fit — i.e. it is
   derived from a fit against **v₀**, converted via H = v₀²/(2g), **not** a fit against
   H = h + d. Expect a small low-speed difference where d is a non-negligible part of H.
   Note also that `step3_fits.m:21` computes `d0mult = b^{3/2}·√g` = **0.6092**, which differs
   from the locked 0.862 by exactly √2. That stored line and the locked constant disagree; the
   locked value was used and `step3_fits.m` was not touched.
6. **d_s could not be reproduced.** `d0_static_20260823_165817.csv` holds 30 rows, all with
   empty `qa_flags`, giving 0.7210 ± 0.1844 — not the locked 0.751 ± 0.170, n = 23. The locked
   value was used on trust; the n = 23 subset is not recoverable from that export. d_s is
   report-only and is not plotted.
7. **4B shows 5 of 54 curves.** Filters are explicit and deterministic (start ≥ −1000 cm/s²,
   then lowest RMSE to k₀z/m; 47 of 54 eligible) and stated on the panel. **All five selected
   curves are low-v₀ (62–118 cm/s)** — an unavoidable consequence of ranking against a straight
   line, since slow drops have the smallest v²/d₁ correction. The subset is a legibility
   choice and is not evidence about how well the fleet collapses; every statistic in this
   report is computed on all 54.
8. The standing limits are unchanged: acceleration is display-grade (windowed derivative),
   only bins ≤ 0.9 cm are usable, intercepts are never k, φ = 0.624 > φ_cps leaves C(z) and
   A(z) degenerate, and geometry is aliased with capture campaign.

**Fig. 4 changes nothing in the scaling, t_stop or ladder figures**, and (k₀, c₀) remain the
manuscript's law parameters.

## 6. Outputs

```
03_RESULTS\_figures\fig_parameter_ladder.pdf / .png
03_RESULTS\_figures\fig_scaling_models.pdf / .png
03_RESULTS\_figures\fig_tstop.pdf / .png
03_RESULTS\_figures\fig_force_decomposition.pdf / .png
```
