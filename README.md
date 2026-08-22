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
streams frames from the video, detects circles, picks the first frame carrying
exactly N markers, tracks from there, and writes `*_tracks.mat` plus QA
overlays. Stage A
deliberately contains **no** kinematics — no smoothing, no velocity, no event
detection — so that a smoothing or model problem can never cost a full batch.

**Frames are not kept (`keepFrames`, default `'none'`).** Stage A decodes each
frame from the video, filters it, detects on it in memory, and discards it. No
`01_FRAMES` folder is written.

```matlab
process_trial('batch', root, outRoot, struct('keepFrames','all'))   % write PNGs as before
```

Exported frames were only ever a disposable cache — `make_video` self-exports to
`tempdir`, `diag_impact_frame` falls back to the raw clip, and Stage B
re-exports the few frames around the impact for QA. The Stage A *record* is the
detections, tracks and scalars, which are unchanged.

The detections are unchanged too, not merely equivalent: the old path wrote a
filtered PNG and read it straight back, and PNG is lossless, so detecting on the
in-memory filtered frame is bit-identical. Both paths call the same
`src/detect_circles_frame.m`.

**Disk implications.** A 12 s clip at ~2.8 kfps is ~34k frames; the autoWindow
subset was already only ~10% of that, and now even those PNGs are gone. What
remains per trial is the detections `.mat`, the tracks, the scalars, and (from
Stage B) a few dozen event-frame PNGs. `01_FRAMES` stops growing entirely
unless you ask for it. Existing `01_FRAMES` folders are still read, but only
under `reuse`/`resume`; `retry` and `overwrite` re-derive from the raw video.

Because the frames are gone, every trial is stamped with `meta.provenance` —
filter, detection parameters, which pass produced them, MATLAB version and the
pipeline git commit — so any frame can be regenerated from raw video + code +
parameters. See [docs/QUANTITIES.md](docs/QUANTITIES.md).

**Raw video location.** The tools that need the raw clip (`diag_impact_frame`,
`make_video`, Stage B's event-frame export) resolve it through
`src/find_raw_video.m`, rooted at `$JERBOA_RAW_ROOT` if set, else
`D:\ME_GRANULAB\Test Batches`. The search is recursive, so the variable need
only name an ancestor of the clips.

**`autoWindow` (on by default).** Campaign-2 clips run ~12 s (~40k frames)
around a ~2k-frame event, so detecting every frame spends nearly all of its time
on empty bed. Before processing, `process_trial` pre-scans the video once with
`VideoReader` — no `imfindcircles` — flagging frames that contain any red pixel
at all (`any(R > 150 & G < 100, 'all')`, the same crude mask `diag_raw_clip`
uses). Only this range is then decoded and detected:

```
[firstRed − pad(1),  lastRed + pad(2)]      clamped to [1, nFrames]
```

`opts.windowPad` defaults to `[200 500]`: 200 frames of pre-impact context for
the velocity-peak search, 500 frames of settling tail for `find_stop` to
extrapolate across. Each trial prints its decision:

```
autoWindow: frames 18240-22310 of 40116 (10.1%)
```

If no frame is flagged, the **full** range is used and `NO_RED_CONTENT` is
warned — a clip is never silently narrowed away. The same fallback applies if
the pre-scan throws.

```matlab
process_trial('batch', root, outRoot, struct('autoWindow', false))          % whole video
process_trial('batch', root, outRoot, struct('windowPad', [400 800]))       % wider margins
```

`autoWindow` narrows the **detection range only**. It never affects which
trials are selected, and it is independent of `dryRun` and `policy`.

**Frame indices are window-relative.** `firstValidFrame` and every tracking
index are relative to the processed window, not the video. The window's absolute
start is recorded as `meta.windowStart` (with `windowEnd` and `autoWindow`) and
in both scalars CSVs, so any index resolves back:

```
absolute video frame = windowStart + firstValidFrame + trackingIndex − 2
```

For a full-range window `windowStart` is 1 and this reduces to the familiar
`firstValidFrame + index − 1`. When a `01_FRAMES` cache *is* reused, its window
is taken from the PNG filenames on disk — `export_frames` names each file by its
absolute index — so a cache written under a different window cannot shift the
offset.

Stage B needs no offset: `kd_kinematics` re-zeroes time at impact and works
entirely inside the tracked array.

**Stage B — kinematics.** `track_tracers_2` drives `kd_kinematics`, which
locates impact and stop and produces depth, velocity and the stopping
acceleration. It follows Katsuragi & Durian (2007): velocity-first, with
`a_stop` read from the fitted acceleration discontinuity.

It also writes the per-trial QA that Stage A no longer leaves behind:
`opts.saveEventFrames` (default true) re-exports `[impact − 20, stop + 20]`
frames from the raw clip into the usual `01_FRAMES` mirror with the usual
raw-indexed naming, and `opts.impactCheck` renders the impact-QA PNG. Both are
skipped under `dryRun`/`preview`, are idempotent, and never fail a trial.

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
GLITCH) plus the manual exclusions in `src/get_manual_exclusions.m`.

That trial count predates the 2026-08-21 integrity review below, which added 20
exclusions and 2 relabels. The live figure is the one `load_kinematics_set`
prints on every load — read it there rather than from this table.

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

## Data notes

### Release-height offsets in the Tight and Wide campaigns

The Tight campaign's release-height reference sat **~+10–13 mm high**: median
`v0 / sqrt(2gh)` = 1.14, 1.07, 1.05 at h = 45, 65, 85 mm, matching a
constant-offset model `sqrt((h + dh)/h)`. The ratio falls toward 1 as h grows,
which is the signature of a fixed offset rather than a scale error — a
mis-calibrated fps or a wrong `mmPerPx` would bias every height by the same
factor.

Height labels are **nominal grouping keys**. All quantitative analysis uses
measured `v0`, which is unaffected by where the release reference sat, so the
offset does not propagate into any fit. The v0 guard's Tight warnings at low
heights reflect this documented offset, **not bad trials**.

The Wide campaign shows a smaller **~3%** offset (ratios 1.03–1.04), also
documented, also immaterial to `v0`-based analysis.

Concretely: the Tight 45 mm trials `T03`, `T05`, `T06`, `T09`, `T10` were
flagged by the free-fall check and deliberately **kept**. They are one
population with their clean 45 mm siblings under this offset, and excluding
them would have removed half a height bin to fix a reference error that no
reported quantity depends on.

### 2026-08-21 integrity review

A matcher was run over every trial, scoring each against every height bin on
`v0`, `d_final` and `t_stop`. It flagged **27 trials** whose recorded height
was inconsistent with their own kinematics. The review resolved them as:

| outcome | n | where it lives |
|---|---|---|
| relabeled | 2 | `src/get_relabel_map.m` |
| unflagged — documented Tight offset, not a defect | 5 | kept; see above |
| excluded | 20 | `src/get_manual_exclusions.m`, blocks A–C |

The two relabels are Default `25mm_T02` and `25mm_T03` → **65 mm**, each a
quantitative match to the 65 mm bin (`d2` = 1.2 and 5.7, margin to the
next-best bin ≥ 19) and confirmed by trace overlay. Relabeling is applied at
**read time only** — no stored file changes, and the trialTag still reads
`25mm_...` so provenance back to the capture session survives.

The matcher script `scripts/relabel_candidates.m` and the CSV it writes,
`relabel_candidates.csv`, are the audit trail for all 27 decisions.

> **`scripts/relabel_candidates.m` is not yet committed.** The maintainer's
> copy is the authoritative one and should be added here as-is. It has
> deliberately not been reconstructed: a re-derived matcher would produce its
> own scores, which is precisely what an audit trail must not do.

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

## Workspace hygiene

Four guards against the failure mode where the tree on disk stops matching what
the analysis thinks it is loading.

**Batch labels are normalized.** A bare number becomes `Batch <n>`. The label is
a directory level, so `5` and `Batch 5` would build two sibling trees for the
same batch and split its trials across both — invisible until an analysis
silently loads half the data. The rewrite is printed when it fires; anything
that is not purely digits is left exactly as typed. Applied to the struct/
headless path and to all three dialogs (batch, single, rerun).

**Duplicate `trialTag`s are a hard error.** `load_kinematics_set` refuses to
build a table where one trial contributes two rows — usually a stale
`_kin_scalars.csv` left under an old batch or model folder. It lists every
duplicated tag with the full path of each source file:

```
2 trialTag(s) appear more than once. One trial must produce exactly one row.

  165mm_T03_full_default  (2 copies)
      ...\03_RESULTS\GB\5\Default Model\...\165mm_T03_full_default_kin_scalars.csv
      ...\03_RESULTS\GB\Batch 5\Default Model\...\165mm_T03_full_default_kin_scalars.csv
```

It stops rather than filtering because a duplicate double-weights that trial in
every per-height mean and every fit, and nothing in the results reveals that it
happened. Which copy is current is a decision for you, not the loader.

**Every committed trial gets an impact-QA image.** Stage B runs
`diag_impact_frame` after each `_kin.mat` is written (`opts.impactCheck`,
default true), producing one PNG per trial in

```
<Root>\03_RESULTS\_batch_logs\impact_checks\
```

so a whole batch can be checked by eye afterwards. The call is wrapped: a QA
failure logs one `WARN` line and never fails a trial whose kinematics are
already on disk. Skipped under `dryRun` and `preview`. Set
`struct('impactCheck', false)` to turn it off.

`diag_impact_frame` itself takes `'Show'` (default true); `'Show',false`
requires `'Save',true` and renders invisibly, which is what the batch path uses.

**Archiving finished work.** `scripts/archive_workspace.m` moves processed
`GB\<batch>\<model>` subtrees out of the workspace and zips them.

```matlab
archive_workspace                      % dry run: prints the full plan
archive_workspace('Confirm', true)     % actually move, zip, verify, delete
```

Dry run is the default and prints every source folder, its file count and size
in GB, the destination, and the zip that would be written. It touches only
`01_FRAMES`, `02_SAVED_DETECTIONS` and `03_RESULTS`; raw videos and existing
`_ARCHIVE` contents are never read for anything but collision checks.

Order matters. Before anything moves, every **Default** `*_kin_scalars.csv` is
copied flat into `_ARCHIVE\reference_default_scalars\` and left
**uncompressed** — that is the baseline for the old-vs-new `impactDistPx`
(−360) comparison, so it has to stay readable without unzipping. Subtrees then
move to `_ARCHIVE\<MODEL>\<BATCH>\<frames|detections|results>\`, one zip is
written per `<MODEL>\<BATCH>`, and the uncompressed copy is deleted **only
after** the zip's entry count is read back and matches the number of files
moved. On any mismatch both copies are kept and a warning is raised — a zip
that cannot be shown to be complete never becomes the only copy. Each run
appends to `_ARCHIVE\ARCHIVE_MANIFEST.txt`.

---

## Rules

- **Do not rename** raw videos, `*_tracks.mat`, folders, or trialTags.
- Files whose names contain a space (Dropbox/OneDrive conflicted copies) are
  ignored by `.gitignore` and must never be committed — MATLAB cannot dispatch
  to them, so they are silently dead weight.
