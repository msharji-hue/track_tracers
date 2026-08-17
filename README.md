# track_tracers

MATLAB codebase for tracking rod/foot markers in high-speed video of granular
impact experiments, and for turning those tracks into penetration kinematics.

The repository is organised on one rule:

> **`src/` is called by other code. `scripts/` is called by a human.**

---

## Pipeline

```
Stage A   tracking        process_trial  ->  process_one_trial
Stage B   kinematics      track_tracers_2  ->  kd_kinematics
Audit                     audit_all_trials
Analysis                  depth_scaling
```

**Stage A — tracking.** `process_trial` walks the video tree, resolves each
clip's metadata, and hands one video at a time to `process_one_trial`, which
exports frames, detects circles, picks the first frame carrying exactly N
markers, tracks from there, and writes `*_tracks.mat` plus QA overlays. Stage A
deliberately contains **no** kinematics — no smoothing, no velocity, no event
detection — so that a smoothing or model problem can never cost a full batch.

**Stage B — kinematics.** `track_tracers_2` drives `kd_kinematics`, which
locates impact and stop and produces depth, velocity and the stopping
acceleration. It follows Katsuragi & Durian (2007): velocity-first, with
`a_stop` read from the fitted acceleration discontinuity.

**Audit.** `audit_all_trials` is the review gate. It builds the per-trial table,
applies the automatic rules, and supports click-selection for manual exclusions:

```matlab
T = audit_all_trials(root,'model','Default','condition','GB/full', ...
                     'clickSelect',true,'write',false,'figures','show');
```

**Analysis.** `depth_scaling` is the single analysis entry point, reading through
`load_kinematics_set`.

```matlab
R = depth_scaling(root,'form','ambroso','condition','GB/full');
```

`form` selects the law under test: `ambroso` (`d = (d0^2 H)^(1/3)`), `powerlaw`
(`d ~ v0^alpha`, alpha fitted), `literature` (published forms head to head), or
`velocity` (`d` vs measured `v0`).

---

## Calibration

Settled values, in `src/get_calibration.m`:

| constant | value |
|---|---|
| `mmPerPx` | `0.1079` (pxPerMm 9.27; 2 mm marker = 18.5 px) |
| `bedX` | `4` — `bedPoint1 = [4,0]`, `bedPoint2 = [4,36]` |
| `g_cm_s2` | `980` |
| `trackTolerancePx` | `25` |

`mmPerPx` is set by direct measurement of the marker diameter. An earlier
approach that inferred it from a pre-impact `g_eff` fit was abandoned: the
camera field of view sits on the bed, so only 0–8 pre-impact frames are
captured — too few to fit the fall parabola. Rail friction is therefore
unmeasured, and impact velocity is taken from the imaging (`v0_cm_s`).

`impactDistPx = -360` px, **standardized across all models and containers** for
the 2026-08 unified campaign.

It is the signed rod-to-bed distance at contact, but it is **not** the impact
frame — it centres the impact *search window*. `kd_kinematics` places a
geometric anchor at `min(abs(rodBedDist_px - impactDistPx))`, then searches
`[anchor − 0.5·wMax, anchor + 2·wMax]` for the peak of the smoothed velocity;
impact is that peak. The trigger only has to bracket contact so a spurious
pre-release maximum cannot win, which is why one value serves every model.

```matlab
calib = get_calibration();                    % bedX 4, trigger -360
calib = get_calibration([], '', 'Tight');     % -360; model recorded only
calib = get_calibration([], '', 'Wide');      % -360; model recorded only
```

Retired: per-model `-370` / `-376.001` / `-409`, which tracked each 3D-printed
foot's rod-to-toe offset — differences well inside the search window — and an
earlier per-container `-290` for `dense`. The `model` and `container` arguments
are still accepted and recorded for provenance, but neither changes the trigger.
Trials whose top marker never reaches −360 may report `ANCHOR_OOR`; the
velocity-peak refinement still locates impact in that case.

**Re-validate if the camera framing changes.** The window is defined in pixels
from the bed line, so a change of lens, working distance or crop can move the
true impact outside it. Run `scripts/diag_impact_frame.m` on a few trials per
model after any framing change.

`impactDistPx` is defined *relative* to whichever bed line is in use, so a
bed-line shift moves the anchor but cancels out of depth, which is measured
relative to the impact frame (`z - z(impact_index)`).

**Frame rate.** `resolve_fps` accepts only `1000–6000` Hz. The lower bound
rejects the 600 fps ffmpeg AVI-muxer artefact; the upper bound covers the
three-model runs, which reach 5200 fps and were silently skipped by the old
5000 ceiling.

**Removed.** The pre-impact effective-gravity and rail-friction estimates are
gone from `kd_kinematics`: the pre-impact window is 0–8 frames, far too short to
fit a free-fall parabola against. `a_stop_cm_s2` replaced them.

---

## Dataset state

Default geometry is validated at **269 trials**:

| condition | trials |
|---|---|
| GB / full | 154 |
| GB / shallow | 27 |
| CHIN / dense | 48 |
| CHIN / as_poured | 40 |

After the automatic rules (NaN depth; `d > 4.00 cm`; GB/shallow `d > 2.50 cm`;
GLITCH) plus 26 manual exclusions.

Substrate properties (`src/get_substrate_properties.m`), finalised as the mean
of five independent preparations:

| condition | φ |
|---|---|
| GB / full | 0.624 ± 0.004 |
| GB / shallow | 0.643 ± 0.009 |
| CHIN / as_poured | 0.280 ± 0.004 |
| CHIN / dense | 0.402 ± 0.004 |

The ± on φ is **random only**; particle densities are treated as exact, so
absolute φ carries an additional systematic that cancels between conditions but
not against literature. Retired single-measurement values: 0.629 / 0.636 /
0.276 / 0.409.

---

## Known open issues

**Campaign 1 Tight and Wide recordings: diagnosed, archived.** *Resolved — kept
here because the archived data must not be re-analysed.*

The clips were not short recordings. **Frames were dropped at capture.** The
external HDD could not sustain the raw stream at 4.4–5.1 kfps, so the camera
wrote only a fraction of the frames it reported. The survey across all 502 clips
found a median inter-frame gap of **~42 camera frames**. What looked like a
17–21 ms clip was a ~1.8 s event sampled at roughly 1/42 of the nominal rate.

The consequence is that **the time base is invalid**, not merely coarse. Gaps are
irregular, so no constant `dt` describes the recording and the frame index cannot
be converted to time. Every downstream quantity — velocity, `v0`, `t_stop`,
`a_stop` — is therefore unrecoverable from these files. Re-processing cannot fix
this; the information was never written to disk.

That dataset is **archived, not analysable**. The survey CSV written by
`scripts/survey_capture_integrity.m` is the record of the diagnosis.

Re-acquired **2026-08** at ~2.7–2.8 kfps to internal SSD, under the capture
protocol below. The model drivers (`run_models`, `run_new_models`, `run_pass2`,
`preflight_new_models`, `review_model_trials`, `clear_model_kinematics`) serve
the new campaign.

**GB/shallow campaign retired 2026-08** pending re-acquisition with verified
labels; archived by maintainer.

---

## Capture protocol

The rules that prevent a repeat of the dropped-frame failure. All three are
mandatory.

1. **Record to the internal SSD only.** Never capture directly to an external
   HDD or a network location. The campaign-1 loss was a sustained-write
   bottleneck, and it is silent — the camera reports the nominal rate whether or
   not the frames reach disk. Move clips to external storage *after* capture.
2. **fps ≤ 2800.** The rig's real rate is ~2700–3100, which resolves the impact
   comfortably. Higher rates buy no useful resolution and reintroduce the write
   bottleneck.
3. **Acceptance-test the first drop of every session** with
   `scripts/diag_raw_clip.m` before recording the rest. It checks the clip's
   actual inter-frame spacing against the nominal rate. If it fails, stop and fix
   the capture path — do not record a session on the assumption that it will be
   repairable later. It will not be.

---

## Naming

**From campaign 2 onward every trial carries a model suffix — including
Default.** Trial tags end `_default`, `_tight` or `_wide`:

```
165mm_T03_full_default      165mm_T03_full_tight      165mm_T03_full_wide
```

Campaign 1 tags for Default carried no suffix (`165mm_T03_full`), which is why
`load_kinematics_set` treats an unsuffixed tag as Default rather than unknown.
Both conventions therefore resolve to the same model, and mixed-campaign loads
work without renaming anything.

Model is part of the analysis group key: `depth_scaling` bins by
`(model, condition, dropHeight_mm)`, so the three geometries are never
averaged together.

---

## Rules

- **Do not rename** raw videos, `*_tracks.mat`, folders, or trialTags.
- Files whose names contain a space (Dropbox/OneDrive conflicted copies) are
  ignored by `.gitignore` and must never be committed — MATLAB cannot dispatch
  to them, so they are silently dead weight.
