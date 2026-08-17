# Reported quantities and their code provenance

Every quantity this pipeline reports, defined exactly as the code computes it,
with the file and function that produces it. Where a definition involves a
choice (a window, a threshold, a protocol convention), the choice is stated.

Sign convention throughout: **positive depth = into the bed**, set by
`rod_displacement` (`src/rod_displacement.m`, `'signConvention'`, default `+1`).

---

## Position: `z`, `depthRod_cm`

`src/rod_displacement.m` → `rod_displacement`, called from
`src/kd_kinematics.m` step 1.

The primary position is the **rigid-rod displacement along the bed normal**,
averaged over all markers visible in each frame:

```matlab
z_rod = rod_displacement(trackedX, trackedY, calib.lineA, calib.lineB, ...
                         calib.lineC, calib.mmPerPx);
```

Each marker's displacement is measured **from its own reference position**, then
averaged — not the mean of raw positions. This matters because markers drop out
mid-impact: a raw positional mean changes value whenever the averaged *set*
changes, injecting a step of order the marker spacing (~34 px) into `z(t)`.
Differentiated twice, that step becomes a delta-function artefact in `a(t)`.
Subtracting each marker's own reference removes it entirely. Averaging over N
visible markers also cuts per-marker detection noise by ~√N.

The reference frame is the first frame where **all** markers are finite,
falling back to the frame with the most visible.

---

## Frame indices are window-relative

`scripts/process_trial.m` → `auto_window` / `apply_auto_window`;
`src/process_one_trial.m` step 2.

Stage A exports a **window** of the video, not necessarily all of it. With
`opts.autoWindow` (default true) the window is the red-marker span padded by
`opts.windowPad` (default `[200 500]`), found by a pre-scan that flags frames
where `any(R > 150 & G < 100, 'all')`. With `autoWindow` off, the window is the
whole video.

Every index downstream — the `detect_circles_per_frame` output, `firstValidFrame`,
and therefore `impact_index` and `stopFrame` — is **relative to that window**.
The window's absolute start is recorded as `meta.windowStart` (alongside
`windowEnd` and `autoWindow`) and in both scalars CSVs, so any stored index can
be resolved back to the video:

```
absolute video frame = windowStart + firstValidFrame + trackingIndex − 2
```

For a full-range export `windowStart` is 1 and this reduces to
`firstValidFrame + trackingIndex − 1`.

This offset matters only where an index is mapped **back to the raw video** —
`src/save_tracks.m` (the `origF` column), `scripts/make_video.m`, and
`scripts/diag_impact_frame.m`. Everything that stays inside the exported set
(`make_annotated_frames`, `make_track_qa`) indexes the PNG list on both sides
and needs no offset.

**No timing quantity depends on it.** `kd_kinematics` builds `t = (0:nF-1)*dt`
over the tracked array and re-zeroes it at impact, so `t_s`, `v0`, `t_stop`,
`a_stop` and `d_final` are unchanged by where the window sits.

---

## `impact_index`

`src/kd_kinematics.m` → main body, step 3.

Two stages: a geometric anchor, then a velocity-peak refinement. **The velocity
peak is the answer**; the anchor only bounds the search.

**1. Geometric anchor.** The signed bed-normal distance of every marker is

```matlab
bedNorm  = sqrt(calib.lineA^2 + calib.lineB^2);
d_px_all = (calib.lineA.*trackedX + calib.lineB.*trackedY + calib.lineC)./bedNorm;
```

The **reference marker** is the one furthest from the bed line on the rod side
(`min(col)`) at the first all-finite frame — chosen by geometry alone. It is
deliberately *not* chosen as "nearest `impactDistPx`": that would place the
marker at the trigger value by construction, collapsing impact onto the
first-detection frame. `rodBedDist_px = d_px_all(refMarkerID,:)` is that
marker's track, reported as `kin.rodBedDist_px`, and `kin.refMarkerID` is its
index.

```matlab
[~, geomImpact] = min(abs(rodBedDist_px - calib.impactDistPx));
```

**2. Search window and velocity peak.**

```matlab
vLite = movmean_omitnan(v_raw, 2*wMin+1);
srch  = max(1,geomImpact-round(0.5*wMax)) : min(nF,geomImpact+round(2*wMax));
[~, li] = max(vLite(srch));
impact_index = srch(1) + li - 1;
```

Impact is the peak of the smoothed velocity: free fall accelerates `v` up to
`v0`, drag then decelerates. The window is asymmetric — half a max-window before
the anchor, two max-windows after — and exists so a spurious pre-release maximum
cannot win.

`wMin`/`wMax` come from `'minWindowMs'` (0.5) and `'maxWindowMs'` (4.0)
converted per trial via `dt`, so the window is a fixed *duration*, not a fixed
frame count, and adapts across the fps spread.

**`impactDistPx` is a search-window centre, not the impact frame.** It is a
standardized `-360` px for all models and containers
(`src/get_calibration.m`), applies to the 2026-08 unified campaign, and must be
re-validated with `scripts/diag_impact_frame.m` if the camera framing changes.
`kd_kinematics` warns (`kd_kinematics:triggerOutOfRange`) when it falls outside
the reference marker's observed range; the velocity peak still locates impact in
that case, and the trial is flagged `ANCHOR_OOR` by `scripts/audit_all_trials.m`.

---

## `t_stop` and `a_stop`

`src/kd_kinematics.m` → local function `find_stop`, following
Katsuragi & Durian (2007), *Nat. Phys.* **3**, 420, Methods.

Three steps:

**1. Zero-crossing.** First frame after impact where `v <= 0`. If none exists the
last frame is used and `kd_kinematics:noStop` is warned.

**2. Rebound exclusion.** The projectile can bounce, so the segment used for the
fit must exclude rebound:

```matlab
reb = max(0, -min(v(cross:end),[],'omitnan'));   % rebound magnitude
thr = rebFactor*reb;
```

Walking back from the crossing, frames are taken only while `v > thr`.
`rebFactor` defaults to **2** (`'rebFactor'` name-value).

**3. Linear extrapolation.** A straight line is fitted to that segment and
extrapolated to `v = 0`:

```matlab
[m, ~]  = fit_slope(t(seg), v(seg));
b       = mean(v(seg)) - m*mean(t(seg));
a_stop  = m;            % KD's acceleration discontinuity
t_stop  = -b/m;
```

So `v(t) = a_stop*(t - t_stop)` near the stop. **`t_stop` is the extrapolated
zero, not the zero-crossing frame** — the two differ, and
`scripts/diag_impact_frame.m` draws both so they can be told apart.
`stopFrame = min(abs(t - t_stop))`, clamped to `[impact_index+1, nF]`.

Reported as `kin.t_stop_s`, `kin.a_stop_cm_s2`, `kin.stopFrame`.

If the segment has fewer than 2 usable frames, `t_stop` falls back to the
crossing time and `a_stop` stays `NaN`.

---

## Velocity: `v`, `v_smooth`, and `v0`

`src/kd_kinematics.m` → steps 2 and 4, local function `adaptive_slope`.

`v_raw` is one central difference of `z_rod` (`local_central_diff`) — KD's
"raw v". `v_smooth` and `a` come from **order-1 adaptive line-segment fits** of
`v_raw`: `v_smooth` is the fitted value at the point (the segment-average
velocity KD plot as open circles) and `a` is the segment slope. Window half-width
grows from `wMin` to `wMax` until the slope standard error meets
`'accelTolRel'`/`'accelTolAbs'` (both 0.005), and is recorded per frame as
`kin.winHalfFrames`.

Windows are clamped to `[impact_index, stopFrame]` so no fit straddles the stop
discontinuity.

```matlab
v0_cm_s = v_smooth(impact_index);
```

**`v0 = 0` for zero-drop trials, by protocol.** `src/load_kinematics_set.m` sets

```matlab
K.v0_meas_cm_s          = K.v0_cm_s;
K.v0_cm_s(K.isZeroDrop) = 0;
```

A `d0` trial is released *from contact*: there is no fall and no impact speed.
What `kd_kinematics` reports for those trials is the velocity peak of a release
transient, so it is kept as `v0_meas_cm_s` for QA — a large value flags a trial
that was dropped rather than released — and never fitted against.

A free-fall cross-check `v0_ff = sqrt(2*g*h_cm)` and the ratio
`v0_ratio = v0_cm_s / v0_ff` are computed in `scripts/depth_scaling.m`. The ratio
is **reported only**; `sqrt(2gh)` is never the fitted predictor, because the
carriage runs on a rail whose friction is unmeasured, so release height is a
control setting rather than a measured input.

---

## `d_final`

`src/kd_kinematics.m` → step 5.

```matlab
depthRod_cm = z_rod - z_rod(impact_index);   % depth zeroed at impact
t_s         = t    - t(impact_index);        % time re-zeroed at impact
d_final_cm  = depthRod_cm(stopFrame);
```

Maximum penetration: the depth at the stop frame, with depth measured relative
to the impact frame.

**Zeroing at impact is what makes depth insensitive to the bed line.** A
bed-line shift moves the geometric anchor (the trigger moves with the line), but
it cancels out of `z - z(impact_index)`. Only the *anchor* is affected, not the
depth. Noted in `src/get_calibration.m`.

**Masking.** After computation, `a` and `a_plus_g` are set `NaN` outside
`[impact_index, stopFrame]` — acceleration is meaningless off that interval — and
`v_smooth`/`depthRod_cm` are cut at `stopFrame + postCapFrames`
(`'postCapMs'`, default from `calib.postCapMs`).

`a_plus_g = -a - g` is the net resistive acceleration in the project convention.

---

## Per-height binning

`scripts/depth_scaling.m` → section 3, "per-height means".

**Group key: `(model, condition, dropHeight_mm)`**, the height rounded to
2 dp in cm (`round(S.h_cm,2)`).

- `model` is in the key because the three feet are different geometries;
  averaging a Tight and a Default trial from the same height would blend two
  projectiles into one point. Fits and figures use the same grouping, so models
  are never pooled back together downstream.
- Grouping is by **drop setting**, which is a *label* for which trials are
  repeats. It never enters a fit. Grouping by rounded `v0` instead splits
  repeats across bins and starves conditions with few settings.

`MinRep` (default **3**) is the minimum trials per group; groups below it are
dropped, and a `(model, condition)` pair with fewer than 3 surviving groups is
skipped with a printed note rather than fitted.

**Representative x and its uncertainty:**

| role | column | definition |
|---|---|---|
| x | `v0_mean` | `mean(g.v0_cm_s)` — mean *measured* impact speed in the group |
| x-uncertainty | `v0_sd` | `std(g.v0_cm_s)` |
| y | `d_mean` | `mean(g.d_final_cm)` |
| y-uncertainty | `d_sd` | `std(g.d_final_cm)` |

The summary also prints the measured `v0` span per `(model, condition)`,
excluding zero-drop trials. Fitted exponents should be read against it: a narrow
span leaves a power law poorly determined however good its R² looks.

Fits are taken on these means, which is what the literature plots (Seguin
average about ten experiments per point). Individual-trial fits are reported as
a robustness column; their R² is limited by real trial-to-trial granular
scatter, not by the model, so the within-height relative SD is reported
alongside.

---

## `dNorm`

`scripts/depth_scaling.m` → section 2.

```matlab
K.dNorm = K.d_final_cm .* sqrt(K.rho_g);
```

Uehara et al. (2003), reproduced as Eq. (2) of Newhall & Durian (2003). For a
single foot the projectile density and diameter are constant, so

$$d\sqrt{\rho_g} = C\,v_0^{2/3},\qquad C \sim 1/\mu$$

`rho_g` is the **bulk** density of the medium (see φ below), not the particle
density. Fitted as `ueh_n` / `ueh_C` in the `literature` form.

---

## `d0`

Two distinct quantities share the name. They are not interchangeable.

**Fitted characteristic length** — `scripts/depth_scaling.m`, section 2:

```matlab
K.d0_cm = sqrt(2*G) * K.d_final_cm.^1.5 ./ K.v0_cm_s;
```

This is Ambroso et al. (2005) `d = (d0²H)^(1/3)` inverted, `d0 = sqrt(d³/H)`,
recast in measured `v0` by substituting `H = v0²/2g`.

*Caveat, stated plainly:* in the original form `d0` is exactly the penetration at
zero drop height, because the `+d` inside `H = h + d` survives as `v0 → 0`.
Dropping to pure `v0` removes that term, so here `d0` is a characteristic length
fitted from `d` and `v0` rather than the h = 0 depth itself. The two agree
closely when penetration is small next to the fall. Reported as `d0_cm`
(`ambroso` form) and `F1_d0_cm` (`literature` form).

**Measured zero-drop depth** — `src/load_kinematics_set.m`, verbose report:
`mean(D0.d_final_cm)` over `isZeroDrop` trials. This is Ambroso's `d+` measured
directly rather than inverted from the drop trials. **Compare the two**; they are
independent estimates of the same length.

A third, unrelated `d0` appears in the `velocity` and `literature` forms as the
**intercept** of `d = a·v0^(2/3) + d0` and `d = d0 + α·v0` respectively. The
latter permits a negative intercept (de Bruyn & Walsh model it as a Bingham
fluid with a yield stress).

---

## φ and `rho_g`

`src/get_substrate_properties.m` — single source of truth, keyed on
`(material, condition)`. Read by `scripts/depth_scaling.m` section 2 and
`src/process_one_trial.m`.

| condition | φ | ρ_particle (g/cm³) | ρ_bulk (g/cm³) |
|---|---|---|---|
| GB / full | 0.624 ± 0.004 | 2.50 | 1.560 |
| GB / shallow | 0.643 ± 0.009 | 2.50 | 1.607 |
| CHIN / as_poured | 0.280 ± 0.004 | 2.35 | 0.658 |
| CHIN / dense | 0.402 ± 0.004 | 2.35 | 0.945 |

- `rho_particle` — true (skeletal) particle density; the φ reference only.
  Soda-lime glass 2.50; Hess pumice / amorphous aluminium silicate 2.35 from SDS
  specific gravity.
- `rho_bulk` — measured bed bulk density (bed mass / bed volume). **This** is
  what enters buoyancy-type terms and `dNorm`; `depth_scaling` uses
  `sub.rho_bulk_g_cm3` as `rho_g`.
- `phi = rho_bulk / rho_particle`.

φ is the mean of **five independent preparations** (pour → level → weigh). The
± is **random only** — particle densities are treated as exact, so absolute φ
carries an additional systematic that cancels between conditions but not against
literature.

Retired single-measurement values: 0.629 / 0.636 / 0.276 / 0.409.

Unknown `(material, condition)` returns `ok = false` with NaNs so callers flag
rather than silently mis-scale; `depth_scaling` errors on it.

---

## Drop height

`dropHeight_mm`, read from `_kin_scalars.csv` by `src/load_kinematics_set.m`.

The value recorded on disk is the physical drop height, so it is used directly
as the height axis with no correction step. `trialTag` is an identifier only.

`scripts/depth_scaling.m` derives `h_cm = dropHeight_mm / 10` and the free-fall
cross-check `v0_ff = sqrt(2*g*h_cm)` from it; see the v0 section above for why
that cross-check is reported but never fitted against.

---

## Frame rate

`src/resolve_fps.m`. Cascade, first plausible value wins:

1. `<tag>_scalars.csv` `fps` — written before `repair_fps` and never rewritten,
   so it holds the authoritative track-time rate
2. `meta.fps_true` — the repaired value
3. `tracks.fps`

Plausible band `LO = 1000`, `HI = 6000` Hz. The lower bound rejects the 600 fps
ffmpeg AVI-muxer artefact; the upper covers the three-model runs. Returns `NaN`
with `src = 'none'` when nothing is plausible, so callers skip the trial rather
than crash on an absurd `dt`.

---

## Removed quantities

`g_eff` and `friction_over_m` were removed from `kd_kinematics`. The camera field
of view sits on the bed, so only 0–8 pre-impact frames are captured — far too few
to fit the curvature of a fall parabola. Rail friction is therefore unmeasured,
and impact velocity is taken from the imaging (`v0_cm_s`). `a_stop_cm_s2`
replaced them.

`mmPerPx` is set by direct measurement of the marker diameter (2 mm = 18.5 px,
pxPerMm 9.27), not inferred from a `g_eff` fit.
