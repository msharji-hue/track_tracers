function make_video(trialTag, style, root, varargin)
% MAKE_VIDEO  Single entry point for trial video generation.
%
%   make_video(trialTag)                    % default style 'tracer'
%   make_video(trialTag,'annotated')        % overlay detections + kinematics
%   make_video(trialTag,'tracer',root)
%
%   style : 'tracer'    -> src/make_tracer_video    (marker trails)
%           'annotated' -> src/make_annotated_video (detections + readout)
%
%   Replaces make_display_video / make_display_video_trial / make_trial_video,
%   which differed only in how they located the trial and which renderer they
%   called. Nothing is picked by hand: the trial is found by tag, the frame
%   range comes from the saved kinematics, and the bed line comes from
%   get_calibration(). The renderers themselves are unchanged.
%
%   OPTIONS (name-value)
%       'RawRoot'    raw video root (default D:\ME_GRANULAB\Test Batches)
%       'PadBefore'  frames before impact (default 40)
%       'PadAfter'   frames after stop    (default 60)
%       'Save'       write an mp4 (default false = preview only)
%       'OutputFps'  output frame rate (default 30)
%       'SlowFactor' repeat each frame this many times (default 2)
%
%   BED LINE. get_calibration returns bedPoint1/bedPoint2 in ORIGINAL
%   (unrotated) pixel coords -- the same frame both renderers expect, and the
%   same frame detect_circles_per_frame reports centres in. The points are
%   passed straight through with no conversion.
%
%   FRAME RANGE. kin.impact_index and kin.stopFrame index the TRACKED array.
%   Tracking frame 1 is firstValidFrame within the EXPORTED WINDOW, and that
%   window starts at meta.windowStart in the raw video, so
%       raw frame = windowStart + firstValidFrame + index - 2
%   (windowStart is 1, and this reduces to the old form, for a full export)
%   Frames are exported starting at the padded impact frame, so the renderers
%   see a 1-based array of their own; 'annotated' therefore receives a kin
%   struct whose event indices have been shifted into EXPORTED-frame
%   coordinates. Trials that were excluded keep their tracks but not their
%   kinematics; for those the full tracked span is used ('tracer' only, since
%   'annotated' has nothing to plot).

if nargin < 2 || isempty(style), style = 'tracer'; end
if nargin < 3 || isempty(root),  root  = 'D:\ME_GRANULAB\JerboaImpact'; end

opt.RawRoot    = 'D:\ME_GRANULAB\Test Batches';
opt.PadBefore  = 40;
opt.PadAfter   = 60;
opt.Save       = false;
opt.OutputFps  = 30;
opt.SlowFactor = 2;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

style    = lower(char(style));
trialTag = char(trialTag);
if ~ismember(style, {'tracer','annotated'})
    error('make_video:badStyle', ...
          'style must be ''tracer'' or ''annotated'', got "%s".', style);
end
fprintf('\n=== %s clip: %s ===\n', style, trialTag);

% ── 1) Locate the trial ───────────────────────────────────────────────────
T = dir(fullfile(root,'03_RESULTS','**',[trialTag '_tracks.mat']));
K = dir(fullfile(root,'03_RESULTS','**',[trialTag '_kin.mat']));
if isempty(T)
    error('make_video:noTracks','No %s_tracks.mat under %s', trialTag, root);
end
tracksPath = fullfile(T(1).folder, T(1).name);
S    = load(tracksPath, 'meta', 'tracks');
meta = S.meta;
fvf  = meta.firstValidFrame;

% Stage A may have exported only a WINDOW of the video (opts.autoWindow), so
% firstValidFrame is an index into that window, not into the video. The
% absolute video frame for tracking index k is
%     windowStart + firstValidFrame + k - 2
% which reduces to the old firstValidFrame + k - 1 when the window started at
% frame 1, i.e. for every full-range export.
wStart = 1;
if isfield(meta,'windowStart') && isfinite(meta.windowStart)
    wStart = meta.windowStart;
end
base = wStart + fvf - 1;          % absolute frame of tracking index 1

haveKin = ~isempty(K);
if haveKin
    Kn  = load(fullfile(K(1).folder, K(1).name), 'kin');
    kin = Kn.kin;
    startFrame = max(1, base + kin.impact_index - 1 - opt.PadBefore);
    endFrame   =        base + kin.stopFrame    - 1 + opt.PadAfter;
    fprintf('impact idx %d, stop idx %d, firstValidFrame %d, windowStart %d\n', ...
            kin.impact_index, kin.stopFrame, fvf, wStart);
else
    if strcmp(style,'annotated')
        error('make_video:noKinematics', ...
              ['%s has no _kin.mat (excluded trial), so the annotated style ' ...
               'has no z(t)/v(t) to plot. Use style ''tracer''.'], trialTag);
    end
    startFrame = base;
    endFrame   = base + meta.nTracked - 1;
    fprintf('no _kin.mat (excluded trial) -- using the full tracked span\n');
end

% ── 2) Resolve fps ────────────────────────────────────────────────────────
% resolve_fps rejects the 600 fps ffmpeg muxer artefact, which would otherwise
% corrupt every time caption and the plotted time axis.
tracksVar = [];
if isfield(S,'tracks'), tracksVar = S.tracks; end
[fps_true, fpsSrc] = resolve_fps(tracksPath, meta, tracksVar);
if ~isfinite(fps_true)
    error('make_video:noFps', ...
          'No plausible fps for %s (resolve_fps returned none).', trialTag);
end
fprintf('fps_true = %.0f (source: %s)\n', fps_true, fpsSrc);

% ── 3) Locate the raw video ───────────────────────────────────────────────
% Matching must be EXACT on the filename stem and must include the material
% and container folders. A substring test on the height label is not enough:
% '25mm' is contained in '125mm' and '325mm', which is how an earlier version
% opened CHIN/as_poured/125mm_T04.avi when asked for 25mm_T04_dense.
V = dir(fullfile(opt.RawRoot,'**','*.avi'));
V = V(~[V.isdir]);
paths = unique(string(fullfile({V.folder}', {V.name}')));  % Windows lists *.avi/*.AVI twice

stem = sprintf('%s_T%02d', meta.heightLabel, meta.trialNum);   % e.g. 25mm_T04
[~, stems] = arrayfun(@(p) fileparts(p), paths, 'UniformOutput', false);
stems = string(stems);

matDir  = string(filesep) + string(meta.material)  + string(filesep);  % \CHIN\
contDir = string(filesep) + string(meta.container) + string(filesep);  % \dense\

hit = strcmpi(stems, stem) & contains(paths, matDir,  'IgnoreCase', true) ...
                           & contains(paths, contDir, 'IgnoreCase', true);
if ~any(hit)
    error('make_video:noMatch', ...
        'No raw video with stem "%s" under a %s%s path in %s', ...
        stem, matDir, contDir, opt.RawRoot);
end

% Prefer the original capture over the transcoded copy: the transcoded AVIs
% carry the 600 fps muxer artefact.
cand = paths(hit);
orig = cand(~contains(cand, "transcoded", 'IgnoreCase', true));
if ~isempty(orig)
    cand = orig;
else
    fprintf('only a transcoded copy exists -- fps is taken from resolve_fps\n');
end
if numel(cand) > 1
    fprintf('%d candidates after filtering, using the first:\n', numel(cand));
    fprintf('   %s\n', cand);
end
videoPath = char(cand(1));
fprintf('video: %s\n', videoPath);

v        = open_video(videoPath);
nAvail   = floor(v.Duration * v.FrameRate);
endFrame = min(endFrame, nAvail);
fprintf('frames %d - %d\n', startFrame, endFrame);

% ── 4) Export sharpened frames ────────────────────────────────────────────
framesDir = fullfile(tempdir, 'display_frames', trialTag);
if exist(framesDir,'dir'), rmdir(framesDir,'s'); end
export_frames(v, startFrame, endFrame, framesDir, 'sharpen');

% ── 5) Detect circles (same defaults as the main pipeline) ────────────────
dp = struct();
dp.radiusRange       = [9 33];
dp.sensitivity       = 0.85;
dp.edgeThresh        = 0.10;
dp.polarity          = 'bright';
dp.alphaG            = 0.50;
dp.betaB             = 0.50;
dp.doCLAHE           = true;
dp.medianK           = 3;
dp.showPreviewEveryN = 0;
dp.showCenters       = true;
dp.heightLabel       = trialTag;
fprintf('Detecting circles...\n');
detectOut = detect_circles_per_frame(framesDir, dp);

% ── 6) Bed line straight from the pipeline calibration ────────────────────
% The impact trigger is per-MODEL, so pass the model through when the trial
% records one; container is recorded for provenance only.
model = '';
if isfield(meta,'model'), model = meta.model; end
container = '';
if isfield(meta,'container'), container = meta.container; end
calib = get_calibration([], container, model);
fprintf('calib: bedX %g, impactDistPx %g, mmPerPx %g\n', ...
        calib.bedX, calib.impactDistPx, calib.mmPerPx);

% ── 7) Shared render params ───────────────────────────────────────────────
rp = struct();
rp.outputFps   = opt.OutputFps;
rp.slowFactor  = opt.SlowFactor;
rp.fps_true    = fps_true;
rp.mmPerPx     = calib.mmPerPx;
rp.bedColor    = [1.00 0.90 0.20];
rp.circLW      = 2;
rp.trailLen    = 8;
rp.trailAlpha  = 0.5;
rp.trailSize   = 6;
rp.rotateCCW   = true;
rp.bedPoint1   = calib.bedPoint1;   % original (unrotated) px, as expected
rp.bedPoint2   = calib.bedPoint2;

fprintf('fps_true=%.0f (file reports %.0f) | outputFps=%d | slowFactor=%d | slowdown=%.0fx\n', ...
    rp.fps_true, v.FrameRate, opt.OutputFps, opt.SlowFactor, ...
    rp.fps_true/opt.OutputFps*opt.SlowFactor);

outDir = fullfile(root,'03_RESULTS','_batch_logs','display_clips');
if opt.Save && ~isfolder(outDir), mkdir(outDir); end

% ── 8) Dispatch to the chosen renderer ────────────────────────────────────
switch style
    case 'tracer'
        rp.markerColor = [0.35 0.95 0.25];
        fprintf('Rendering preview...\n');
        make_tracer_video(framesDir, detectOut, rp, false);
        if opt.Save
            rp.outputName = fullfile(outDir, [trialTag '_tracer.mp4']);
            make_tracer_video(framesDir, detectOut, rp, true);
        end

    case 'annotated'
        rp.markerColor = [0.20 0.80 1.00];   % cyan
        rp.circLW      = 3;
        rp.trailAlpha  = 0.6;

        % Renderer contract: det.detect.{centersCell,radiiCell}, indexed over
        % the EXPORTED frames.
        det = struct('detect', detectOut);

        % Shift the event indices from tracked-array coords into exported-frame
        % coords. Tracked index i is raw frame (base + i - 1); exported frame 1
        % is raw frame startFrame, so exported = base + i - startFrame.
        nExported     = numel(dir(fullfile(framesDir,'*.png')));
        kinShift      = kin;
        kinShift.impact_index = base + kin.impact_index - startFrame;
        kinShift.stopFrame    = base + kin.stopFrame    - startFrame;
        rp.startFrame = 1;
        rp.stopFrame  = nExported;
        fprintf('exported %d frames | impact -> %d, stop -> %d (exported coords)\n', ...
                nExported, kinShift.impact_index, kinShift.stopFrame);

        fprintf('Rendering preview...\n');
        make_annotated_video(framesDir, det, kinShift, rp, false);
        if opt.Save
            rp.outputName = fullfile(outDir, [trialTag '_annotated.mp4']);
            make_annotated_video(framesDir, det, kinShift, rp, true);
        end
end

if ~opt.Save
    fprintf('Preview only. Re-run with ''Save'',true to write the mp4.\n');
end
end
