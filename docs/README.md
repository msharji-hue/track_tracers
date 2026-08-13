# track_tracers

MATLAB codebase for detecting and tracking tracer markers (rod/foot markers) in high-speed videos from granular impact experiments. The goal is stable marker IDs across frames, robust handling of missed detections, and standardized outputs for depth, velocity, acceleration, and downstream modeling/plots.

## Quickstart

The **official entry point** is `scripts/run_track_tracers.m`. 
Never call `track_tracers()` directly from the command window without a `params` struct.
```matlab
% In MATLAB, from the repo root:
run('scripts/run_track_tracers.m')

% Or, with custom params:
params = struct();
params.material    = "GB";
params.batchName   = "Batch1";
params.heightLabel = "H06";
params.trialNum    = 1;
params.version     = "short";
results = track_tracers(params);
```

**Video file resolution:** The tracker accepts `.mov`, `.mp4`, or `.avi` and will
pattern-match `<heightLabel>_<trialID>_<version>.*` if the exact filename is missing.
If no video is found, the error message prints the directory searched, the pattern used,
and the first 30 filenames in that folder to help diagnose path issues.

## Expected workflow
1. Provide a video clip (short test clip recommended).
2. Run `track_tracers.m` (or the main runner script).
3. Inspect `debug/` overlays if tracking looks off.
4. Export tracked trajectories for post-processing.

## Notes
- Bed/surface reference can be applied using a fixed line model when the camera setup is standardized across experiments.
- Keep video paths out of the repo if they are large; store clips locally and reference them in a config file.
