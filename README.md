# track_tracers

MATLAB codebase for detecting and tracking tracer markers (rod/foot markers) in high-speed videos from granular impact experiments. The goal is stable marker IDs across frames, robust handling of missed detections, and standardized outputs for depth, velocity, acceleration, and downstream modeling/plots.

## What’s inside
- Marker detection (thresholding/segmentation + centroid extraction)
- Frame-to-frame association (ID stability, gating, re-acquisition)
- Debug mode (overlay frames + logs) to diagnose silent tracking failures
- Standardized outputs (pixel coordinates and optional conversion to physical units)

## Expected workflow
1. Provide a video clip (short test clip recommended).
2. Run `track_tracers.m` (or the main runner script).
3. Inspect `debug/` overlays if tracking looks off.
4. Export tracked trajectories for post-processing.

## Notes
- Bed/surface reference can be applied using a fixed line model when the camera setup is standardized across experiments.
- Keep video paths out of the repo if they are large; store clips locally and reference them in a config file.
