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

`impactDistPx` is the signed rod-to-bed distance at contact and depends on the
foot **model**, because each 3D-printed foot has a different rod-to-toe offset:

| model | `impactDistPx` |
|---|---|
| Default | `-370` |
| Tight | `-376.001` |
| Wide | `-409` |

```matlab
calib = get_calibration();                    % Default: bedX 4, trigger -370
calib = get_calibration([], '', 'Tight');     % -376.001
calib = get_calibration([], '', 'Wide');      % -409
```

The `container` argument is accepted and recorded for provenance but **does not**
change the trigger. An earlier per-container value of −290 for `dense` was
tested and retired in favour of one uniform trigger. Dense trials may therefore
report `ANCHOR_OOR`, since the top marker often does not reach −370 in those
trials; the velocity-peak refinement still locates impact in that case.

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

## GB/shallow height reversal

**GB/shallow drop-height labels are recorded in reverse order on disk.**

Nothing is renamed. Filenames, folder names, trialTags and `meta.dropHeight_mm`
all keep the original label. `src/true_drop_height.m` is the single source of
truth for the correction:

| labelled | actual | | labelled | actual |
|---|---|---|---|---|
| 25 mm | 365 mm | | 285 mm | 125 mm |
| 65 mm | 325 mm | | 325 mm | 65 mm |
| 125 mm | 285 mm | | 365 mm | 25 mm |
| 165 mm | 165 mm (unchanged) | | | |

Only GB/shallow is altered; every other condition passes through unchanged.

The map is a swap and therefore self-inverse, so it must be applied **exactly
once**, at read time, from `true_drop_height` only. Never hardcode it elsewhere.

- Use `dropHeight_true_mm` for **all** physics — v0 checks, scaling, every plot axis.
- Use `dropHeight_mm` and `trialTag` **only** as file identifiers.

`load_kinematics_set` applies the correction for you and returns
`dropHeight_true_mm` alongside a `heightCorrected` flag.

**Evidence.** With the labelled heights, GB/shallow `v0 / sqrt(2gh)` reached
3.9–4.2 — impossible, since a ratio above 1 means faster than free fall — and
measured v0 was nearly constant (~275, 265, 240 cm/s) across a 5× range of
labelled heights. After correction every ratio falls to ~1.0–1.1.

Depth (`d_final_cm`) is **not** affected by this correction.

Markers on disk: `03_RESULTS/HEIGHT_LABELS_README.txt` and the per-folder
`_HEIGHT_REVERSED.txt` files, both written by `scripts/write_height_labels.m`.

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

**Tight and Wide GB/full recordings are unusable.** Those clips are ~17–21 ms
long against Default's ~1.8 s, so the impact event is not contained in them.
Those datasets are not analysable until re-acquired or reprocessed. The model
drivers (`run_models`, `run_new_models`, `run_pass2`, `preflight_new_models`,
`review_model_trials`, `clear_model_kinematics`) are kept pending that
investigation.

---

## Rules

- **Do not rename** raw videos, `*_tracks.mat`, folders, or trialTags.
- Do not hardcode the height map anywhere outside `true_drop_height`.
- Files whose names contain a space (Dropbox/OneDrive conflicted copies) are
  ignored by `.gitignore` and must never be committed — MATLAB cannot dispatch
  to them, so they are silently dead weight.
