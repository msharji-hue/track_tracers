function [z_cm, nUsed, refFrame, perMarker_cm] = rod_displacement(trackedX, trackedY, ...
                                                  lineA, lineB, lineC, mmPerPx, varargin)
% ROD_DISPLACEMENT  Rigid-rod displacement along the bed-normal axis, averaged
% over ALL markers visible in each frame — composition-invariant.
%
%   The markers sit on the RIGID ROD attached to the foot (none of them is the
%   foot/toe itself). Because the rod is rigid, every marker undergoes the SAME
%   displacement, so we average each marker's displacement FROM ITS OWN
%   REFERENCE POSITION rather than averaging raw positions.
%
%   Why this matters: a raw positional mean, mean(pos(visible,t)), changes value
%   whenever a marker appears/disappears, because the SET being averaged changes.
%   In these trials markers drop out mid-impact (e.g. frames 45/56/70), which
%   would inject a STEP of order the marker spacing (~34 px) directly into z(t).
%   Differentiating a step twice yields a delta-function artefact in a(t) — i.e.
%   exactly the unphysical (a+g) spikes we are trying to avoid. Subtracting each
%   marker's own reference first removes this entirely: dropping a marker changes
%   which terms are averaged, but not the average's value.
%
%   Averaging over N visible markers also reduces per-marker detection noise
%   by ~sqrt(N), and adapts naturally to however many markers survive in a given
%   trial (no fixed subset, no per-trial tuning).
%
%   [z_cm, nUsed, refFrame, perMarker_cm] = rod_displacement(...)
%
%   Inputs
%     trackedX/Y  [nMarkers x nFrames] tracked positions (NaN where unmatched)
%     lineA/B/C   bed line  A*x + B*y + C = 0   (from get_calibration)
%     mmPerPx     calibration scale
%   Name-value
%     'refFrame'  reference frame index (default: first frame where ALL markers
%                 are finite; falls back to the frame with the most visible)
%     'minMarkers' minimum visible markers for a frame to yield a value (default 1)
%     'signConvention' +1 or -1 (default +1). Depth is returned so that
%                 POSITIVE z = motion toward/into the bed.
%
%   Outputs
%     z_cm         [1 x nFrames] mean rod displacement (cm), NaN where < minMarkers
%     nUsed        [1 x nFrames] how many markers contributed to each frame
%     refFrame     reference frame actually used
%     perMarker_cm [nMarkers x nFrames] each marker's own displacement (cm)

    p = inputParser;
    addParameter(p,'refFrame',[]);
    addParameter(p,'minMarkers',1);
    addParameter(p,'signConvention',1);
    parse(p,varargin{:});
    o = p.Results;

    [nM, nF] = size(trackedX);

    % ── signed perpendicular distance of every marker to the bed line ────
    bedNorm = sqrt(lineA^2 + lineB^2);
    d_px    = (lineA .* trackedX + lineB .* trackedY + lineC) ./ bedNorm;

    % ── reference frame: prefer the first frame with ALL markers visible ──
    if isempty(o.refFrame)
        nVis = sum(isfinite(d_px), 1);
        refFrame = find(nVis == nM, 1, 'first');
        if isempty(refFrame)
            [~, refFrame] = max(nVis);          % fallback: most-visible frame
        end
    else
        refFrame = o.refFrame;
    end

    % ── each marker's displacement from ITS OWN reference ────────────────
    ref = d_px(:, refFrame);                    % [nM x 1]
    perMarker_px = d_px - ref;                  % NaN stays NaN
    % markers with no reference contribute nothing
    perMarker_px(~isfinite(ref), :) = NaN;

    % ── composition-invariant mean over visible markers ──────────────────
    finiteMask = isfinite(perMarker_px);
    nUsed      = sum(finiteMask, 1);
    dispZeroed = perMarker_px;
    dispZeroed(~finiteMask) = 0;                % mask already excludes these
    sumDisp    = sum(dispZeroed, 1);
    z_px       = sumDisp ./ nUsed;
    z_px(nUsed < o.minMarkers) = NaN;

    % ── to cm, with sign so that positive = into the bed ─────────────────
    z_cm         = o.signConvention .* z_px  .* mmPerPx ./ 10;
    perMarker_cm = o.signConvention .* perMarker_px .* mmPerPx ./ 10;
end
