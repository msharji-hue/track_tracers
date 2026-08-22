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

## Frames are derived data

`src/detect_circles_stream.m`; `src/process_one_trial.m` step 2.

**Stage A does not keep frames.** With `keepFrames = 'none'` (the default) each
frame is decoded from the video, filtered, detected on, and discarded. No
`01_FRAMES` folder is written.

Frames are a *cache*, not a record. The Stage A record is the detections, the
tracks and the scalars. Any frame can be reproduced exactly from three things:

| input | where |
|---|---|
| the raw video | the capture tree (`$JERBOA_RAW_ROOT`) |
| the code | `meta.provenance.gitCommit`, `matlabVersion` |
| the parameters | `meta.provenance.filterType`, `detectParams`, `detectPass`, `windowStart`, `windowEnd` |

`meta.provenance` is stamped on every trial:

```
frameSource    stream | png-cache | png-export | detections-cache
keepFrames     none | all
filterType     the filter applied before detection (apply_filter)
detectParams   the detection parameters actually used
detectPass     pass1 | pass2-backup  (which pass produced these detections)
windowStart/End the exported window, absolute video frames
matlabVersion  version()
gitCommit      short sha, '-dirty' if the tree had uncommitted changes
processedOn    timestamp
```

A `-dirty` commit means the sha does not pin what ran; treat those results as
unversioned.

**Detections are unchanged by streaming.** The old path wrote a filtered PNG
and read it back; PNG is lossless and `apply_filter` returns uint8 RGB, so
detecting on the in-memory filtered frame is bit-identical. Both paths call the
same `src/detect_circles_frame.m` per frame, so they cannot drift apart.

**Where frames still exist.** Stage B re-exports a small subset around the
impact (`opts.saveEventFrames`, default `[20 20]` frames either side) into the
usual `01_FRAMES` mirror with the usual raw-indexed naming, so visual QA has
something to read without decoding the clip again. `make_video` self-exports to
`tempdir`, and `diag_impact_frame` falls back to the raw video when no PNG is
present.

---

## Frame indices are window-relative

`scripts/process_trial.m` → `auto_window` / `apply_auto_window`;
`src/process_one_trial.m` step 2.

Stage A processes a **window** of the video, not necessarily all of it. With
`opts.autoWindow` (default true) the window is the red-marker span padded by
`opts.windowPad` (default `[200 500]`), found by a pre-scan that flags frames
where `any(R > 150 & G < 100, 'all')`. With `autoWindow` off, the window is the
whole video.

Every index downstream — the detection output, `firstValidFrame`,
and therefore `impact_index` and `stopFrame` — is **relative to that window**.
The window's absolute start is recorded as `meta.windowStart` (alongside
`windowEnd` and `autoWindow`) and in both scalars CSVs, so any stored index can
be resolved back to the video:

```
absolute video frame = windowStart + firstValidFrame + trackingIndex − 2
```

For a full-range window `windowStart` is 1 and this reduces to
`firstValidFrame + trackingIndex − 1`.

This offset matters only where an index is mapped **back to the raw video** —
`src/save_tracks.m` (the `origF` column), `scripts/make_video.m`, and
`scripts/diag_impact_frame.m`. Everything that stays inside the exported set
(`make_annotated_frames`, `make_track_qa`) indexes the window on both sides and
needs no offset — `make_annotated_frames` maps window index to video frame
internally when it reads from the video rather than from a PNG cache.

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

**Saved field names differ from the internal ones.** `v_smooth` and `z_rod` are
local variables inside `kd_kinematics`; they are *not* struct fields. What
reaches disk is:

| internal | saved as |
|---|---|
| `v_smooth` | `kin.v` |
| `depthRod_cm` | `kin.depthRod_cm`, aliased `kin.z` |
| `a` | `kin.a` |
| `a_plus_g` | `kin.a_plus_g` |

Reading `kin.v_smooth` or `kin.z_smooth` gets you nothing — that was a real bug
in `make_annotated_video`. Consumers want `kin.v` and `kin.z`.

```matlab
v0_cm_s = v_smooth(impact_index);
```

**Zero-drop trials are quarantined.** See the section below; `v0_cm_s` is 0 by
protocol on those rows and the raw value is kept as `v0_meas_cm_s` for QA.

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

### `a + g` — net grain acceleration

`kin.a_plus_g = g - a`, the acceleration the grains supply once gravity is
accounted for. Depth is positive into the bed, so a free-falling rod has
`a = +g` and three checkpoints fix both the sign and the offset:

| state | `a` | `a + g` |
|---|---|---|
| free fall | `+g` | `0` — nothing but gravity |
| at rest | `0` | `+g` — the bed carries the weight |
| decelerating | `-abs(a)` | `g + abs(a)` — large positive |

This matches KD 2007 Fig. 1.

**Corrected 2026-08.** The formula was previously `-a - g`: the same quantity
offset by a constant `-2g`, with the rest state inverted, so free fall sat at
`-2g` and a resting rod at `-g` instead of `+g`.

**It was display-only.** No fit and no reported scalar ever read it — `a_stop`
is fitted to velocity in `find_stop`, `v0` comes from `v_smooth`, `d_final` from
`depthRod_cm`, and neither `load_kinematics_set` nor `depth_scaling` touches it.
Only figures were wrong.

**Stored `_kin.mat` files written before the fix still carry the old formula**
in their `a_plus_g` column. That is harmless: the raw `a` trace beside it was
never affected, and every consumer now recomputes through `src/net_accel.m`
(`g - a`) rather than reading the stored column. Already-processed trials
therefore render correctly with no Stage B re-run. Re-running Stage B refreshes
the column itself, but nothing depends on that.

---

## Zero-drop trials are quarantined

`src/load_kinematics_set.m`, immediately after `isZeroDrop` is set.

A `d0` trial is released **from contact**: no fall, no impact, no deceleration.
`kd_kinematics` is a *dynamic-impact* pipeline and is correct for its domain —
h = 0 is simply outside it. Everything it derives from the impact/stop event is
therefore **undefined** on those rows, not merely noisy:

| column | zero-drop value | why |
|---|---|---|
| `v0_cm_s` | `0` | protocol: no fall, so no impact speed |
| `d_final_cm` | `NaN` | depth at a "stop" located by a zero-crossing that never belonged to a deceleration |
| `t_stop_s` | `NaN` | extrapolated from that same non-event |
| `a_stop_cm_s2` | `NaN` | the slope of a fit to it |

The raw numbers are preserved in `*_meas` companion columns —
`v0_meas_cm_s`, `d_final_meas_cm`, `t_stop_meas_s`, `a_stop_meas_cm_s2` — for
**QA only**. A large `v0_meas_cm_s` flags a trial that was dropped rather than
released. Never fit against them; they are diagnostics, not measurements.

`kd_kinematics` is deliberately **not** changed. Its machinery is right for
impacts; the fix is to stop asking it about a case it does not model.

**Backlog.** `d0` needs a dedicated quasi-static extraction: the settled
penetration under the foot's own weight, measured with no impact event. Until
that exists these rows carry no usable depth, and the fitted `d0` from the
`ambroso` form has no independent measurement to be checked against.

**Nothing downstream can fit them.** `scripts/depth_scaling.m` excludes them
three times over — the `isZeroDrop` test, `v0 > 0`, and `isfinite(d_final_cm)` —
and then *asserts* `~any(K.isZeroDrop)`. The redundancy is the point: if a
future edit relaxes any one condition, the assert fails loudly rather than
quietly fitting a released-from-contact trial as though it were a drop. The
loader's NaN-depth exclusion rule stays exempt for zero-drop rows, or the
quarantine would delete the very rows it is meant to preserve.

---

## Figure 1 — kinematics figure

`scripts/fig_kinematics.m`. Depth, velocity and net acceleration against time
from impact, per foot model.

**Pre-agreed caption language:**

> Curves are pointwise medians over replicate drops. Zero-drop trials are
> reserved for scalar measurements and are not shown in the time-history
> figures. Net-acceleration curves begin at the free-fall value a+g = 0 at the
> impact instant and terminate at the rest value a+g = g (dotted line) at each
> height's unified stop; both endpoints are physically determined boundary
> states (free fall: a = g; rest: a = 0); the steep terminal drop is the
> measured stopping discontinuity (Delta a = a+g_terminal - g). Values between
> are the measured ensemble medians. All quantitative fits use unmodified
> kinematics.

**Caption language for `ColorBy` `'height'` (the default):**

> Colors indicate nominal drop height on the ladder shared by all three
> geometries; measured v0 for each condition appears in the velocity panel.
> Nearby heights differ in impact speed by less than the trial-to-trial spread
> at the upper ladder, so height labels are nominal grouping keys; all
> quantitative analysis uses measured v0.

**Display layer.** Every step below affects the drawing only. Nothing is
written back, and no fit sees any of it:

| step | applies to |
|---|---|
| resample each replicate onto a uniform grid (`GridMs`, default 0.2 ms) relative to impact | all rows |
| after each trial's OWN detected `t_stop` (`kin.t_stop_s`, never the quarantined table scalar), fill ALL of that trial's remaining grid points — to the end of the grid, not just to its own last sample — with its measured **rest state** — `v = 0`, depth = its `d_final`, `a + g = NaN`. Filling to the grid end is what keeps a trial in the pointwise aggregate for the whole drawn interval; stopping at its last sample leaves NaNs *inside* the region and can void the median mid-curve. Where a recording ends before its detected `t_stop`, rest-fill begins at that last sample instead — an approximation, counted and reported as `nShortRec` under `'Diagnose', true` | depth, v (not `a + g`: the pipeline computes no post-stop acceleration, and none is invented here) |
| pointwise `Average` across the **full replicate set** (`'median'` default, or `'mean'`) at every grid point | all rows |
| curve terminates at the height's **median replicate `t_stop`** — the same endpoint in all three panels | all rows |
| `movmean` per row (`SmoothMs`: depth 1.0 ms, v 2.0 ms, a+g 1.5 ms), applied AFTER aggregation; NaN-aware, with the pre-smoothing NaN mask re-applied afterward | each row independently |
| **edge-preserving** smoothing — window half-width = distance to the edge of the finite support, so the first and last samples are returned unchanged and the window stays symmetric | `a + g` only |
| before impact, draw only where **support ≥ `MinReplicates`** — trials carry different amounts of pre-impact recording, so the median steps as the contributing set changes | all rows |
| after impact, draw the computed median only where **support ≥ `MinReplicates`**, cutting the thin-support tail. Not a clamp: nothing is floored at `g` | `a + g` only |
| append a final segment from the last drawn sample to `(t_stop, 0)` — the rest state, closing the small-positive end the median has by construction | `v` only |
| prepend a leading segment from `(t = 0, a + g = 0)` to the first computed sample — the free-fall boundary state | `a + g` only |
| append one final segment from the last computed value down to `a + g = g` at the unified stop — the rest boundary state | `a + g` only |
| on the default linear axis, nothing is masked or clipped. Under `'AccelScale','log'` only, mask `a + g <= AccelYMin` to NaN, since a log axis cannot render non-positive values | `a + g` only |

**Why rest-extension replaced "terminate below `MinReplicates`".** Replicates
stop at different times. Under the old rule, past the earliest stop the
aggregate was taken over only the longest-lasting trials — survivor bias — and
the curve was pulled toward them before ending early. Filling each trial's rest
state keeps every replicate in the aggregate for its whole duration; run with
`'Diagnose', true` for the per-height replicate `t_stop` spread, the support
count at the curve's endpoint under both rules, and confirmation that display
smoothing does not shift the `a + g` onset. `MinReplicates` (default 3) now
gates only whether a height is **plotted at all**, not where a curve ends.

A consequence to check in the render: `v` decays to ~0 at the endpoint by
**construction of the median** — about half the replicates are already at rest
there — not by being forced to zero.

**Zero-drop trials are excluded by default** (`'ShowZeroDrop'`, default
`false`): they are reserved for scalar measurements and are not shown in these
time-history figures. The individual-intrusion overlay rule (a zero-drop trial
drawn as its own thin curve, depth/v only, when its depth range over the whole
recorded trace exceeds `IntrusionThreshCm`) is retained behind the same option
rather than deleted, and is otherwise unchanged: a threshold rule evaluated
against the unwindowed trace, not a hardcoded tag list.

**`a + g` is computed over `[impact, stop]` only, and anchored at both ends.**
`kd_kinematics` masks its stored `a` — and therefore `a + g` — to NaN outside
`[impact_index, stopFrame]`, and the figure plots exactly that between its two
boundary anchors. Nothing inside that interval is reconstructed: no pre-impact
acceleration is estimated from the velocity or from anywhere else.

The two anchors are the only drawn values not taken from the trace, and both
are **exact consequences of the definition** `a + g = g - a`, not estimates:

| boundary | state | `a` | `a + g` |
|---|---|---|---|
| onset, as `t → 0⁻` | free fall | `g` | **0** |
| stop, at `t_stop` | at rest | `0` | **`g`** |

Neither appears in the stored trace, because the pipeline computes `a` only
between impact and stop. That makes them unrecorded, not unknown — the same
status the rest state has always had here. A leading segment runs from
`(0, 0)` to the first computed sample; a trailing segment runs from the last
computed sample down to `g`.

This is a deliberate reversal of the earlier "no synthetic onset point" rule,
and the distinction is what makes it one. That rule rejected inserting a
*chosen* value (`1e1`) to force a visible rise. Zero is not chosen: it is what
`g - a` equals for every falling trial at every height. Anchoring to it means
the steep rise reads as the measured jump from its true starting state, instead
of the curve materialising at `~g` partway up.

**`a + g` (row c) is on a linear y-axis by default.** Over `[impact, stop]` this
rig's values span roughly `1e3`–`2.2e4` cm/s² — about one decade — where a
linear axis reads more directly than a log one. (KD 2007 Fig. 1c is logarithmic
because their 20 µs sampling resolves `1e1`–`1e4`, three decades.) Pass
`'AccelScale','log'` for the KD-matched log presentation, floored at
`AccelYMin` (default `1e2` cm/s²). Depth and v are never clamped; v's y-axis
lower limit is fixed at 0.

The `a + g` y-range is **one shared range across all three geometry figures**:

```
ymin = 0                    hard floor, never raised
ymax = 1.05 × (global max over every plotted a+g curve, all models)
```

The floor is a hard 0 rather than a padded data minimum because 0 is the
free-fall boundary value the onset segments anchor to — padding downward would
waste the axis, raising it would cut them off. Only the top gets headroom. Both
limits are printed in the run summary.

**Axis margins.** Depth and velocity are padded by 5% of their drawn data
range, so no curve, marker or endpoint sits on the frame. Limits are computed
once from the pooled drawn data and shared across models within each row.

**The `a + g` curve ends in the rest state.** One final segment runs from the
last computed value down to `a + g = g` at the unified stop: at rest the bed
carries the rod's weight, so `a = 0` and `a + g = g` exactly. The drop across
the segment is the per-height **acceleration discontinuity** — how much net
upward acceleration the grains were still supplying at the instant the rod
stopped — and a light dotted reference line at `g` is drawn behind the data so
the drops can be seen landing on a marked baseline. The line is unlabelled: it
is there to be seen against, not read off.

It terminates at the **same** `t` as the depth and velocity curves — one stop
per `(model, height)`, shared by that figure's three panels, never a separately
computed acceleration endpoint and never shared between geometries. Nothing is
drawn past it. `'Diagnose', true` asserts the three endpoints are identical,
comparing where each row is *last drawn* — for the two anchored rows that is the
end of their terminal segment, which the `a + g` support gate may place beyond
the end of the computed series.

**Velocity is anchored the same way, to its own rest value.** At the unified
stop the median `v` ends small-**positive** by construction: about half the
replicates are already at rest at 0 there and the rest are still moving, so the
median sits just above zero. A final segment to `(t_stop, 0)` closes that —
same helper, same endpoint. Unlike the `a + g` drop it carries no physical
content; it removes a construction artefact rather than displaying a measured
discontinuity, and it does not assert the rod decelerated faster than measured.

The three panels do **not** share a *start*, and are not meant to: depth and
velocity begin at `-PreCapMs`, `a + g` at its free-fall anchor at `t = 0`. The
assertion is on endpoints only, deliberately.

**Both ends are support-gated.** A row is drawn only where at least
`MinReplicates` replicates contribute — the same threshold that decides whether
a height is plotted at all.

*Before impact this applies to all three rows.* Trials do not carry the same
amount of pre-impact recording: Stage A's window starts where the red markers
appear, so at −10 ms a height may have three trials contributing where at −2 ms
it has ten. The median then steps as the contributing **set** changes, which is
the kink that was visible in the Default velocity curve near −5 ms — an artefact
of who is in the average, not of the motion. Each row uses its own support
count, since the three need not thin out together.

Neither gate can punch an interior hole. Support is non-increasing away from
`t = 0` in both directions — going back, a trial drops out at the start of its
own recording and never returns; going forward, at its own stop — so each gate
trims a contiguous head or tail. `'Diagnose', true` prints both cut times per
height.

**After impact only the `a + g` tail needs gating.** Depth and velocity are rest-extended, so
every replicate contributes at every grid point. `a + g` is not — the pipeline
computes no post-stop acceleration — so past the earliest replicate stop the
pointwise median is taken over fewer and fewer trials, until it is wandering
over two or three. That produced a dip to ~665 cm/s², *below* `g`, just before
a terminal drop: not physics, a median with almost nothing left in it.

The computed curve is therefore drawn only where **support ≥ `MinReplicates`**
— the same threshold that already decides whether a height is plotted at all —
and the terminal segment runs from that last well-supported sample to
`(unified stop, g)`. The unified stop is unchanged, so the endpoint assertion
is unaffected. The gate always removes a contiguous **tail**, never a mid-curve
hole: after impact each trial contributes over `[0, its own stop]`, so the
count is non-increasing in `t`. It is applied to the aggregate *before*
smoothing, so the smoother's trailing edge sits at the gate rather than
reaching past it into samples that will not be drawn.

**This is a gate, not a floor.** Nothing is clamped to `g` anywhere. `a + g`
below `g` mid-penetration is real and is kept exactly as computed: `a + g < g`
is simply `a > 0`, the grain force momentarily below the foot's own weight,
which is ordinary during penetration. Only the thin-support tail is removed.

What *would* be a defect is the median dipping below `g` while most of its own
replicates stay above it — an aggregate driven by one or two outliers rather
than by the population. So `'Diagnose', true` reports, for each height whose
curve dips below `g` after impact, the fraction of the replicates contributing
at the dip minimum whose own `a > 0`, and flags it below **0.70**. It is a
diagnostic, not a filter: it decides nothing about what is drawn.

The `'trials'` QA style applies **no** gate: a single trial has support 1
everywhere, so the ensemble threshold would blank the entire row in the one
view whose purpose is showing what an individual trial recorded.

**Why the `a + g` row alone is smoothed edge-preserving.** `movmean`'s
`'shrink'` rule truncates the window at a boundary but still averages about
half a window there. Depth and velocity begin `PreCapMs` before impact, where
the curve is flat, so that costs nothing. `a + g` begins *at* impact, on the
steepest part of the trace: at the impact sample `v` is maximal, so `a ≈ 0` and
`a + g ≈ g` (~0.1e4). Averaging that measured onset together with the rise above
it reported the mid-rise value (0.3–1.6e4) as the curve's starting point.

Setting the window's half-width to the distance from the edge of the finite
support returns the first and last samples unchanged while keeping the window
symmetric, so nothing is phase-shifted. `'Diagnose', true` asserts the first
`a + g` sample equals the aggregate's.

**Colour encodes one of two things (`'ColorBy'`).**

| mode | index | why |
|---|---|---|
| `'height'` (default) | rank of the nominal drop height on the ladder shared by all three geometries — 25 mm coldest, 365 mm warmest | a **rank on a fixed ladder** makes the three figures step through identical colours *by construction*, whatever `v0` each geometry reached at a given height |
| `'v0'` | median impact speed through one shared range | the earlier behaviour, retained unchanged |

Nothing here modifies a measured `v0`. It remains the plotted data in the
velocity panel, the `medV0` on every curve, and the number every analysis uses.

Zero-drop curves are excluded from both references even when drawn: their `v0`
is 0 by quarantine rather than by measurement, which would anchor the `v0`
scale at zero and put a phantom rung at the foot of the ladder. A height not on
the ladder — a zero-drop overlay at `h = 0` — falls to the nearest rung, taking
the coldest colour.

The ladder is the **union** of the heights across geometries, not the
intersection, so a height present in only one geometry keeps its own rung and
does not shift the colours of the heights around it in the others.

**Both references are built from the full loaded set, before the `'Models'`
filter** — not from the curves that were built. That distinction is the fix for
a real defect in the `'v0'` mode: `fig_kinematics` filters its table by
`'Models'` immediately after loading, so a range taken from the built curves was
a property of the *invocation*, not of the dataset. Drawing all three geometries
in one call gave one mapping; regenerating a single geometry on its own — the
natural way to redraw one figure — gave another. Nothing was wrong inside the
per-model loop; the range handed to it was already wrong.

Both are reproduced without loading any traces: group the full set by
`(model, height)`, take the median `v0` per group, keep groups with at least
`MinReplicates` trials. That group set is a slight *superset* of the curves
actually drawn — a group can clear the threshold here and still lose its curve
if a trial proves unreadable. For `'v0'` that can only widen the range a little;
for `'height'` it can only add a rung no curve occupies, which shifts nothing,
since a height's colour is its rank and the ladder is the same in every figure
either way. Both are stable across invocations, which is the property being
bought.

The active reference is printed in the run summary and returned as `R.clim` /
`R.ladder` (both populated in either mode, alongside `R.colorBy`). Each figure's
extreme curves are printed under `colour check` so the claim is verifiable from
the log: in `'height'` mode the same height must print the same index in every
figure; in `'v0'` mode the same `v0` must, and every line must quote the same
`clim`.

**Style.** A **black** dashed vertical line at `t = 0` in every panel —
impact is the reference the whole figure is drawn against, so it is not light
grey chart furniture; a solid thin line at `depth = 0` in the depth panel (the
bed surface); a light dotted line at `g` in the `a + g` panel. No others.

Axes convention (`hold`/`grid`/`box`) comes from `src/apply_fig_style.m`,
shared with `scripts/depth_scaling.m` so the two figure families stay matched;
change the house style there and both follow.

**Ticks.** Major ticks are set explicitly and every one is labelled. Time uses
the smallest round step from 5, 10, 20, 25, 50 ms that does not crowd the axis
— a `grid3x3` panel is a third the width of a per-model one and gets a
correspondingly lower tick budget. The `y` axes use the standard round step
(`1`, `2`, `2.5`, `5` × 10ⁿ). No minor ticks: `apply_fig_style` leaves them at
MATLAB's default off, and enabling them here would desynchronise this figure
from `depth_scaling`. Every panel keeps full numeric labels on both axes and
its own x-label.

**Pre-impact window.** `'PreCapMs'`, default **10 ms**. Keep it an exact
multiple of `GridMs` (10/0.2 = 50) so `t = 0` lands *on* a grid point — that is
what makes the first ensemble `a + g` sample the impact sample itself rather
than a point already partway up the rise.

**Layout.** `'per-model'` (default) writes three figures, one per model, each
a 3×1 vertical stack sharing the time axis, single-column APS width. `'grid3x3'`
reproduces the earlier combined rows-×-models figure in one file. Axis limits
are identical across models within each row in both layouts.

**`Style`** selects the presentation: `'mean'` (default) is the figure above;
`'trials'` draws every kept trial individually and unaveraged, for QA (no rest
extension — each trial ends at its own stop), sharing the same drawing code so
the two are directly comparable.

**Reportable scalars.** `fig1_kinematics_stops.csv` (model, height, median v0,
median `t_stop`, n) is written alongside the figures for the cross-geometry
comparison.

---

## Which trials are excluded, and by what

`src/load_kinematics_set.m` → the rules block; `src/get_manual_exclusions.m`.

Two kinds of exclusion, kept deliberately apart.

**Automatic rules** — conditions the code can test, applied to `h > 0` trials
only (a released-from-rest trial has no fall or deceleration to test):

| flag | rule | option |
|---|---|---|
| `NaNDEPTH` | `d_final_cm` not finite | `'dropNaNDepth'` (default true) |
| `GLITCH` | \|v\| does not fall monotonically from `v0` to zero between impact and stop; an increase above 15% of `v0` marks a tracking artefact | `'dropGlitch'` (default true) |
| `depth cut` | `d_final_cm` greater than the cut | `'depthCut'` (default off) |

**Manual exclusions** — judgements the code cannot make, listed in
`src/get_manual_exclusions.m`, one tag per line with its own reason and date:

```matlab
"345mm_T02_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
```

That file is the single source of truth. `load_kinematics_set` applies it **by
default**, so every analysis reaching data through the loader drops the same
trials — this is what stopped the per-script copies from drifting apart when
`load_default_gb` was retired.

The list is kept as one **block per review decision**, concatenated at the end
of the function. The block header carries the full reasoning; each line still
carries a short reason and a date, so a tag stays self-describing if a block is
ever split. Two reviews are recorded:

| block | n | reason |
|---|---|---|
| 2026-08-20 | 23 | stop-fit failure or a bad velocity profile (all Default GB/full) |
| 2026-08-21 A | 5 | `v0`/`t_stop` inconsistent with any ladder height — off-ladder drop or failed trial |
| 2026-08-21 B | 14 | measurement failure: `v0` inconsistent with the trial's own `d_final` and/or above any free-fall bound |
| 2026-08-21 C | 1 | `v0` matches the 125 mm class but `d_final` = 1.71 cm sits outside that bin's core; no clean match |

Block A's five are not individually odd. As a *group* they share a coherent
low-`v0` / long-`t_stop` signature matching no rung of the ladder — and
specifically not 65 mm, where the two relabeled 25 mm trials went. A consistent
signature that fits nothing is a different failure from a noisy trace, and it
is why they are excluded rather than relabeled: there is no bin to move them to.

**Overriding.** Passing `'exclude'` **replaces** the list rather than adding to
it:

```matlab
load_kinematics_set(root)                                       % reviewed list
load_kinematics_set(root, 'exclude', strings(0,1))              % none
load_kinematics_set(root, 'exclude', [get_manual_exclusions(); "extra_tag"])
```

`depth_scaling` forwards `'exclude'` only when the caller actually set it — its
own default is `[]`, meaning *unset*. Forwarding an unset default would have
replaced the reviewed list with an empty one and silently re-included every
excluded trial.

**Traceability.** The returned table carries `keep` and `reason`, so an
exclusion can be traced rather than silently applied, and the verbose report
prints the per-rule counts plus any listed tag that matched no file on disk —
a stale tag left after a rename is otherwise invisible.

Anything expressible as a rule belongs in the rules, not in the manual list.
To retire a manual exclusion, delete its line: git history is the record of
what was once excluded, and a commented-out tag reads as ambiguous forever.

**Excluding is not the only answer.** A trial whose *height label* is wrong but
whose data are sound belongs in `get_relabel_map` (below), not here: excluding
throws away a good measurement, relabeling keeps it in the bin it came from.
Exclusion is for trials with no bin to go to.

---

## Height relabeling: `dropHeight_mm` vs `dropHeight_asRecorded`

`src/get_relabel_map.m`, applied by `src/load_kinematics_set.m`.

A few trials carry a recorded height label that review established is wrong.
`get_relabel_map` returns a `trialTag` → corrected `dropHeight_mm` table, one
row per trial with its reason and date, and the loader applies it **at read
time**.

| column | meaning |
|---|---|
| `dropHeight_mm` | the **corrected** height — what every axis, bin and fit should use |
| `dropHeight_asRecorded` | the label as it came off the CSV; provenance, not a measurement |
| `relabeled` | true where `dropHeight_mm` came from the map |

**Nothing on disk is modified.** Raw videos, `01_FRAMES`, tracks, `_kin.mat`
and `_kin_scalars.csv` all keep their original labels — including the height
baked into the `trialTag` itself, so a relabeled trial reads `25mm_...`
forever. That is deliberate: renaming files would break provenance back to the
capture session and to every diagnostic ever run against that tag. The cost is
that anything reading the CSVs *directly* rather than through
`load_kinematics_set` sees the original label.

Applied **before** `isZeroDrop` and before the zero-drop quarantine, because
the corrected height is the physical truth: were a relabel ever to move a trial
to or from `h = 0`, every rule downstream must see the corrected value.

As of 2026-08-21 the map holds **two** entries — Default `25mm_T02` and
`25mm_T03` → **65 mm** — each a quantitative match to the 65 mm bin on `v0`,
`d_final` and `t_stop` (`d2` = 1.2 and 5.7, margin to the next-best bin ≥ 19),
confirmed by trace overlay. The loader prints a summary line whenever a relabel
fires, and reports any mapped tag that matched no file.

---

## The `v0` consistency guard

`src/load_kinematics_set.m` → the guard block, after exclusions.

Every kept `h > 0` trial is checked against the free-fall speed its label
implies, `sqrt(2gh)` with `g` = 980 cm/s² (as `calib.g_cm_s2`). The measured
ratio is returned as `v0_freefall_ratio`.

| side | test | option | what it means |
|---|---|---|---|
| high | `v0 > v0RatioMax · sqrt(2gh)` | `'v0RatioMax'` (default 1.15) | physically impossible — the **label** or the **time base** is wrong |
| low | `v0 < 0.5 · sqrt(2gh)` | fixed | a spoiled or snagged release (e.g. the known 85 mm Default trial at 39.9 cm/s) |

**It warns and counts; it never drops.** Which trials to remove is a review
decision that belongs in `get_manual_exclusions` with a reason attached, not in
an automatic filter that would quietly reshape a bin the next time a threshold
moved. The two sides are not symmetric: exceeding free fall is impossible once
tolerance is allowed for, while falling short of it merely needs explaining, so
the low-side threshold is deliberately loose — half of free fall is far below
any honest air-drag loss over these heights.

Real drops land a few percent *below* free fall, so the healthy ratio is just
under 1. Two campaigns sit systematically above it because their release-height
reference was offset; those warnings are documented behaviour, not bad trials —
see the data notes in `README.md`.

Option names on `load_kinematics_set` are matched case-insensitively and an
unrecognised name is an **error**. A mistyped `'v0RatioMax'` would otherwise be
absorbed as a struct field nothing reads, disabling the guard with no sign that
it had happened.

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

**Measured zero-drop depth** — *not currently available.* This would be
Ambroso's `d+` measured directly rather than inverted from the drop trials, and
the natural independent check on the fitted `d0`. The zero-drop trials are
quarantined (see below), so there is nothing to compare against until the
quasi-static extraction exists.

A third, unrelated `d0` appears in the `velocity` and `literature` forms as the
**intercept** of `d = a·v0^(2/3) + d0` and `d = d0 + α·v0` respectively. The
latter permits a negative intercept (de Bruyn & Walsh model it as a Bingham
fluid with a yield stress).

---

## `A_bare` and `A_hull` — intruding foot area

`scripts/fig_foot_schematic.m`, computed from the STLs in `cad/`.

Both are areas of the foot **projected normal to the drop axis**, counting only
the part that enters the bed.

**The cut.** The STL drop axis is `−Y` (toes first). Everything at
`y <= CutY` intrudes; `CutY` defaults to **−50 mm**, the top of the rectangular
beam (foot-bar junction). Triangles crossing that plane are clipped
(Sutherland–Hodgman) rather than dropped —
discarding them would remove area exactly at the plane where it is being
measured.

**The projection.** Clipped triangles are projected onto **XZ**, the plane
normal to the drop axis. Projecting onto XY instead would give the side
silhouette, which is a different quantity.

| quantity | definition |
|---|---|
| `A_bare` | area of the union of the projected triangles, as a `polyshape`. Gaps between the toes are **holes** and are not counted. |
| `A_hull` | area of the convex hull of the same projected vertices. The swept envelope, gaps included. |

Both constants are **foot-only**. The inclined bar and the marker post begin
at `y = −50` and enter the bed once the toe tip is deeper than `z = 1.12 cm`;
beyond that, the bare projected area of everything below the surface grows by
about **0.10 cm² per mm** (2.12 → 2.96 cm² at `z = 2.0 cm`, 3.73 cm² at
`z = 2.9 cm`). The depth dependence is written by the script to
`foot_area_vs_depth.csv` and plotted in `fig_foot_area_vs_depth`.

Reported in cm² (`polyshape` area is mm², divided by 100).

**Why both.** The three models differ *only* in toe splay. `A_bare` is therefore
identical across them — same beam, same toes, just rotated — and only `A_hull`
grows. Any difference in penetration between models cannot be attributed to
intruding area in the bare sense; the hull is what changes.

| model | `A_bare` (cm²) | `A_hull` (cm²) | ratio |
|---|---|---|---|
| Tight | 2.122 | 2.607 | 1.23 |
| Default | 2.122 | 3.495 | 1.65 |
| Wide | 2.122 | 4.052 | 1.91 |

Computed independently (trimesh/shapely, 2026-08-20) and asserted by the script
at 2% relative tolerance on every run. A mismatch aborts before any figure is
written.

**Sensitivity to the cut.** Both areas change by roughly **4–5% per mm** of
`CutY`, and the shift is nearly uniform across the three models — so the
**ratios above are insensitive within about ±0.5 mm of the cut**, even though
the absolute values are not. Quote `A_hull/A_bare` when comparing models; quote
an absolute area only alongside the `CutY` it was computed at. Every value the
script writes carries its `CutY`, in the CSV and in the figure title.

They are **not** insensitive to counting the bar: `A_hull/A_bare` =
1.23 / 1.65 / 1.91 at the cut becomes 1.15 / 1.45 / 1.66 at `z = 2.9 cm`.

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
