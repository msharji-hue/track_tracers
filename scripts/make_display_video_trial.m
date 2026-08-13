function make_display_video_trial(trialTag, varargin)
%MAKE_DISPLAY_VIDEO_TRIAL  Tracer clip for ONE named trial, no prompts.
%
%   Same flow and same helpers as make_display_video.m, but nothing is picked
%   by hand: the trial is found by tag, the frame range comes from the saved
%   kinematics, and the bed line comes from get_calibration().
%
%   BED LINE. get_calibration returns bedPoint1/bedPoint2 in ORIGINAL
%   (unrotated) pixel coords -- the same frame make_tracer_video expects, and
%   the same frame detect_circles_per_frame reports centres in. So the calib
%   points are passed straight through with no conversion, exactly as the
%   hand-clicked points were.
%
%   FRAME RANGE. kin.impact_index and kin.stopFrame index the TRACKED array,
%   which starts at meta.firstValidFrame in the raw video:
%       raw frame = firstValidFrame + index - 1
%   Trials that were excluded keep their tracks but not their kinematics; for
%   those the full tracked span is used instead and a note is printed.
%
%   USAGE
%       make_display_video_trial('165mm_T06_dense')
%       make_display_video_trial('165mm_T06_dense','Save',true)
%
%   OPTIONS
%       'Root'       results root   (default D:\ME_GRANULAB\JerboaImpact)
%       'RawRoot'    raw video root (default D:\ME_GRANULAB\Test Batches)
%       'PadBefore'  frames before impact (default 40)
%       'PadAfter'   frames after stop    (default 60)
%       'Save'       write an mp4 (default false = preview only)
%       'OutputFps'  output frame rate (default 30)
%       'SlowFactor' repeat each frame this many times (default 2)

opt.Root       = 'D:\ME_GRANULAB\JerboaImpact';
opt.RawRoot    = 'D:\ME_GRANULAB\Test Batches';
opt.PadBefore  = 40;
opt.PadAfter   = 60;
opt.Save       = false;
opt.OutputFps  = 30;
opt.SlowFactor = 2;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

trialTag = char(trialTag);
fprintf('\n=== display clip: %s ===\n', trialTag);

% ── 1) Locate the trial ───────────────────────────────────────────────────
T = dir(fullfile(opt.Root,'03_RESULTS','**',[trialTag '_tracks.mat']));
K = dir(fullfile(opt.Root,'03_RESULTS','**',[trialTag '_kin.mat']));
if isempty(T)
    error('make_display_video_trial:noTracks', ...
          'No %s_tracks.mat under %s', trialTag, opt.Root);
end
S = load(fullfile(T(1).folder, T(1).name), 'meta');
meta = S.meta;
fvf  = meta.firstValidFrame;

if ~isempty(K)
    Kn = load(fullfile(K(1).folder, K(1).name), 'kin');
    startFrame = max(1, fvf + Kn.kin.impact_index - 1 - opt.PadBefore);
    endFrame   =        fvf + Kn.kin.stopFrame    - 1 + opt.PadAfter;
    fprintf('impact idx %d, stop idx %d, firstValidFrame %d\n', ...
            Kn.kin.impact_index, Kn.kin.stopFrame, fvf);
else
    startFrame = fvf;
    endFrame   = fvf + meta.nTracked - 1;
    fprintf('no _kin.mat (excluded trial) -- using the full tracked span\n');
end

% ── 2) Locate the raw video ───────────────────────────────────────────────
% Matching must be EXACT on the filename stem and must include the material
% and container folders. A substring test on the height label is not enough:
% '25mm' is contained in '125mm' and '325mm', which is how an earlier version
% opened CHIN/as_poured/125mm_T04.avi when asked for 25mm_T04_dense.
V = dir(fullfile(opt.RawRoot,'**','*.avi'));
V = V(~[V.isdir]);
paths = unique(string(fullfile({V.folder}', {V.name}')));   % Windows lists *.avi/*.AVI twice

stem = sprintf('%s_T%02d', meta.heightLabel, meta.trialNum);   % e.g. 25mm_T04
[~, stems] = arrayfun(@(p) fileparts(p), paths, 'UniformOutput', false);
stems = string(stems);

matDir  = string(filesep) + string(meta.material)  + string(filesep);  % \CHIN\
contDir = string(filesep) + string(meta.container) + string(filesep);  % \dense\

hit = strcmpi(stems, stem) & contains(paths, matDir, 'IgnoreCase', true) ...
                           & contains(paths, contDir, 'IgnoreCase', true);
if ~any(hit)
    error('make_display_video_trial:noMatch', ...
        'No raw video with stem "%s" under a %s%s path in %s', ...
        stem, matDir, contDir, opt.RawRoot);
end

% Prefer the original capture over the transcoded copy: the transcoded AVIs
% carry the 600 fps muxer artefact, which would corrupt the time caption.
cand = paths(hit);
orig = cand(~contains(cand, "transcoded", 'IgnoreCase', true));
if ~isempty(orig)
    cand = orig;
elseif numel(cand) >= 1
    fprintf('only a transcoded copy exists -- fps is taken from meta, not the file\n');
end
if numel(cand) > 1
    fprintf('%d candidates after filtering, using the first:\n', numel(cand));
    fprintf('   %s\n', cand);
end
videoPath = char(cand(1));
fprintf('video: %s\n', videoPath);

v = open_video(videoPath);
nAvail   = floor(v.Duration * v.FrameRate);
endFrame = min(endFrame, nAvail);
fprintf('frames %d - %d\n', startFrame, endFrame);

% ── 3) Export sharpened frames ────────────────────────────────────────────
framesDir = fullfile(tempdir, 'display_frames', trialTag);
if exist(framesDir,'dir'), rmdir(framesDir,'s'); end
export_frames(v, startFrame, endFrame, framesDir, 'sharpen');

% ── 4) Detect circles (same defaults as the main pipeline) ────────────────
params = struct();
params.radiusRange       = [9 33];
params.sensitivity       = 0.85;
params.edgeThresh        = 0.10;
params.polarity          = 'bright';
params.alphaG            = 0.50;
params.betaB             = 0.50;
params.doCLAHE           = true;
params.medianK           = 3;
params.showPreviewEveryN = 0;
params.showCenters       = true;
params.heightLabel       = trialTag;
fprintf('Detecting circles...\n');
detectOut = detect_circles_per_frame(framesDir, params);

% ── 5) Bed line straight from the pipeline calibration ────────────────────
calib = get_calibration();
fprintf('calib: bedX %g, impactDistPx %g, mmPerPx %g\n', ...
        calib.bedX, calib.impactDistPx, calib.mmPerPx);

% ── 6) Render params ──────────────────────────────────────────────────────
rp = struct();
rp.outputFps   = opt.OutputFps;
rp.slowFactor  = opt.SlowFactor;
if isfield(meta,'fps_true') && isfinite(meta.fps_true) && meta.fps_true > 1000
    rp.fps_true = meta.fps_true;      % pipeline value; immune to the muxer cap
else
    rp.fps_true = v.FrameRate;
end
rp.markerColor = [0.35 0.95 0.25];
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

% ── 7) Preview, then optionally save ──────────────────────────────────────
fprintf('Rendering preview...\n');
make_tracer_video(framesDir, detectOut, rp, false);

if opt.Save
    outDir = fullfile(opt.Root,'03_RESULTS','_batch_logs','display_clips');
    if ~isfolder(outDir), mkdir(outDir); end
    rp.outputName = fullfile(outDir, [trialTag '_tracer.mp4']);
    make_tracer_video(framesDir, detectOut, rp, true);
else
    fprintf('Preview only. Re-run with ''Save'',true to write the mp4.\n');
end
end
