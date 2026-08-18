function R = diag_impact_frame(trialTag, varargin)
%DIAG_IMPACT_FRAME  Visual verification of impact detection for ONE trial.
%
%   Shows the frame kd_kinematics called impact, with the tracked markers, the
%   reference marker, the bed line and the impact search-window centre drawn on
%   it, alongside the velocity trace that actually located the event. Use it to
%   confirm by eye that the located frame really is contact.
%
%   READ-ONLY. Loads *_tracks.mat and *_kin.mat and displays; it never writes to
%   the results tree. Passing 'Save' writes a PNG to the batch-log folder only.
%
%   THIS IS THE RE-VALIDATION TOOL FOR impactDistPx. That value is a standardized
%   -360 px for every model and only centres the SEARCH WINDOW: kd_kinematics
%   anchors at min(abs(rodBedDist_px - impactDistPx)) and then takes the peak of
%   the smoothed velocity inside [anchor - 0.5*wMax, anchor + 2*wMax]. The window
%   is in pixels from the bed line, so a change of lens, working distance or crop
%   can move the true impact outside it. After any framing change, run this on a
%   few trials per model and check that:
%     - the marker overlay sits at the bed at the displayed frame
%     - the velocity peak is inside the window, not at its edge
%     - rodBedDist_px at impact is not wildly far from impactDistPx
%
%   USAGE
%       diag_impact_frame('165mm_T06_full_default')
%       diag_impact_frame(tag,'Root',root,'Model','Tight')
%       R = diag_impact_frame(tag,'Save',true);
%
%   OPTIONS (name-value)
%       'Root'      results root   (default D:\ME_GRANULAB\JerboaImpact)
%       'RawRoot'   raw video root. Defaults to $JERBOA_RAW_ROOT if set, else
%                   D:\ME_GRANULAB\Test Batches. Searched RECURSIVELY, so it
%                   only has to be an ancestor of the clips, not the exact
%                   batch folder. Set the environment variable once to point
%                   the whole toolchain at the current campaign.
%       'Model'     ''|'Default'|'Tight'|'Wide'. Only used when the tracks file
%                   carries no saved calib; otherwise the SAVED calib wins, so
%                   what is drawn is what the trial was actually processed with.
%       'Pad'       frames either side of [impact, stop] in panel B (default 20)
%       'Save'      write a PNG to <Root>/03_RESULTS/_batch_logs/impact_checks
%                   (default false)
%       'Show'      display the figure (default true). Set false for
%                   unattended batch QA: the figure is built invisibly, written
%                   to the PNG, and closed, so a long batch does not accumulate
%                   hundreds of open windows. Show=false requires Save=true --
%                   otherwise the work would be discarded with nothing to see
%                   and nothing on disk, which is never what was meant.
%
%   FRAME SOURCE. Exported frames are a disposable intermediate and are often
%   purged. This looks for the PNG first, under 01_FRAMES (honouring
%   JERBOA_FRAMES_ROOT), and falls back to reading the frame straight out of the
%   raw .avi. Either way the frame shown is the RAW frame
%       rawFrame = meta.windowStart + meta.firstValidFrame + kin.impact_index - 2
%   because kin.impact_index indexes the TRACKED array, firstValidFrame indexes
%   the EXPORTED WINDOW, and windowStart places that window in the video.
%
%   Returns R with the located indices, the drawn geometry and the frame source,
%   so a caller can batch this without re-deriving any of it.

opt.Root    = 'D:\ME_GRANULAB\JerboaImpact';
opt.RawRoot = local_default_raw_root();
opt.Model   = '';
opt.Pad     = 20;
opt.Save    = false;
opt.Show    = true;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

trialTag = char(trialTag);

fprintf('\n=== diag_impact_frame: %s ===\n', trialTag);
fprintf('  Root    : %s\n', opt.Root);
fprintf('  RawRoot : %s\n', opt.RawRoot);
fprintf('  Pad     : %d frames\n', opt.Pad);
fprintf('  Save    : %s\n', local_tern(opt.Save,'yes','no (display only)'));
fprintf('  Show    : %s\n', local_tern(opt.Show,'yes','no (headless QA)'));
if ~opt.Show && ~opt.Save
    error('diag_impact_frame:nothingToDo', ...
        ['Show=false requires Save=true. With neither, the figure would be ' ...
         'rendered and discarded: nothing displayed and nothing written.']);
end

% ── 1) Locate and load ───────────────────────────────────────────────────
T = dir(fullfile(opt.Root,'03_RESULTS','**',[trialTag '_tracks.mat']));
if isempty(T)
    error('diag_impact_frame:noTracks','No %s_tracks.mat under %s', trialTag, opt.Root);
end
tracksPath = fullfile(T(1).folder, T(1).name);
S = load(tracksPath);
meta   = S.meta;
tracks = S.tracks;

K = dir(fullfile(opt.Root,'03_RESULTS','**',[trialTag '_kin.mat']));
if isempty(K)
    error('diag_impact_frame:noKinematics', ...
        ['No %s_kin.mat under %s. This trial has tracks but no kinematics ' ...
         '(excluded, or Stage B not run), so there is no impact frame to check.'], ...
        trialTag, opt.Root);
end
Kn  = load(fullfile(K(1).folder, K(1).name), 'kin');
kin = Kn.kin;

% Calibration: the SAVED one wins, so the drawn geometry is what the trial was
% actually processed with, not what today's get_calibration would return.
if isfield(S,'calib') && isstruct(S.calib) && isfield(S.calib,'impactDistPx')
    calib = S.calib;  calibSrc = 'tracks file';
else
    container = ''; if isfield(meta,'container'), container = meta.container; end
    model = opt.Model;
    if isempty(model) && isfield(meta,'model'), model = meta.model; end
    calib = get_calibration([], container, model);
    calibSrc = 'get_calibration (no calib saved in tracks)';
end

impact = kin.impact_index;
stopF  = kin.stopFrame;
refID  = kin.refMarkerID;
rbAtImpact = kin.rodBedDist_px(impact);

fprintf('\n  calib source     : %s\n', calibSrc);
fprintf('  model / container: %s / %s\n', ...
    local_get(calib,'model','(none)'), local_get(calib,'container','(none)'));
fprintf('  bed line         : x = %g   (A=%.4g B=%.4g C=%.4g)\n', ...
    local_get(calib,'bedX',NaN), calib.lineA, calib.lineB, calib.lineC);
fprintf('  impactDistPx     : %g px  (search-window centre)\n', calib.impactDistPx);
fprintf('  rodBedDist_px    : %.2f px at impact  -> offset %+.2f px from centre\n', ...
    rbAtImpact, rbAtImpact - calib.impactDistPx);
fprintf('  impact_index     : %d   stopFrame: %d   refMarkerID: %d\n', ...
    impact, stopF, refID);
fprintf('  v0_cm_s          : %.2f      d_final_cm: %.4f      t_stop_s: %.5f\n', ...
    kin.v0_cm_s, kin.d_final_cm, kin.t_stop_s);

rbFin = kin.rodBedDist_px(isfinite(kin.rodBedDist_px));
if ~isempty(rbFin)
    inRange = calib.impactDistPx >= min(rbFin) && calib.impactDistPx <= max(rbFin);
    fprintf('  reference marker observed range [%.1f, %.1f] px -> centre %s\n', ...
        min(rbFin), max(rbFin), local_tern(inRange,'REACHABLE','OUT OF RANGE (ANCHOR_OOR)'));
end

% ── 2) The frame ─────────────────────────────────────────────────────────
fvf = meta.firstValidFrame;
% Stage A may have exported only a window of the video, so firstValidFrame is
% an index into that window. windowStart is 1 for a full-range export.
wStart = 1;
if isfield(meta,'windowStart') && isfinite(meta.windowStart)
    wStart = meta.windowStart;
end
rawFrame = wStart + fvf + impact - 2;
[img, frameSrc] = local_get_frame(tracksPath, meta, trialTag, rawFrame, opt);
fprintf('  frame            : raw #%d (windowStart %d + firstValidFrame %d + impact %d - 2)\n', ...
    rawFrame, wStart, fvf, impact);
fprintf('  frame source     : %s\n', frameSrc);

% ── 3) Geometry to draw ──────────────────────────────────────────────────
% Signed distance is d = (A*x + B*y + C)/norm, so the locus at distance d is the
% bed line translated by d along the unit normal (A,B)/norm.
n2   = sqrt(calib.lineA^2 + calib.lineB^2);
uhat = [calib.lineA, calib.lineB] / n2;
bp1  = calib.bedPoint1(:).';
bp2  = calib.bedPoint2(:).';
% Extend the drawn segments across the frame so they are visible whatever the
% recorded endpoints were; the endpoints only record provenance.
dirv = (bp2 - bp1); if norm(dirv) == 0, dirv = [0 1]; end
dirv = dirv / norm(dirv);
L    = max(size(img,1), size(img,2));
bedA = bp1 - 2*L*dirv;   bedB = bp1 + 2*L*dirv;
trgA = bedA + calib.impactDistPx*uhat;
trgB = bedB + calib.impactDistPx*uhat;

mx = tracks.trackedX(:,impact);
my = tracks.trackedY(:,impact);

% ── 4) Figure ────────────────────────────────────────────────────────────
if opt.Show
    fig = figure('Color','w','Position',[60 60 1280 560]);
else
    fig = figure('Color','w','Position',[60 60 1280 560],'Visible','off');
end
tl  = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
title(tl, sprintf('%s  --  impact frame check', strrep(trialTag,'_','\_')), ...
      'FontWeight','bold');

% Panel A: the frame
ax1 = nexttile(tl);
imshow(img,'Parent',ax1); hold(ax1,'on');
plot(ax1,[bedA(1) bedB(1)],[bedA(2) bedB(2)],'-', ...
     'Color',[1.00 0.90 0.20],'LineWidth',1.8,'DisplayName','bed line');
plot(ax1,[trgA(1) trgB(1)],[trgA(2) trgB(2)],'--', ...
     'Color',[0.20 0.80 1.00],'LineWidth',1.6, ...
     'DisplayName',sprintf('search centre (%g px)',calib.impactDistPx));
plot(ax1,mx,my,'o','MarkerSize',8,'MarkerEdgeColor',[0.35 0.95 0.25], ...
     'LineWidth',1.5,'LineStyle','none','DisplayName','tracked markers');
if isfinite(refID) && refID >= 1 && refID <= numel(mx)
    plot(ax1,mx(refID),my(refID),'s','MarkerSize',15, ...
        'MarkerEdgeColor',[1 0.2 0.2],'LineWidth',2.2,'LineStyle','none', ...
        'DisplayName',sprintf('reference marker (#%d)',refID));
end
legend(ax1,'Location','southoutside','Box','off','NumColumns',2);
title(ax1, sprintf('raw frame %d  (tracked index %d)', rawFrame, impact));

txt = sprintf(['impact\\_index = %d\nrodBedDist\\_px = %.1f\n' ...
               'impactDistPx = %g\noffset = %+.1f px\nv0 = %.1f cm/s'], ...
               impact, rbAtImpact, calib.impactDistPx, ...
               rbAtImpact - calib.impactDistPx, kin.v0_cm_s);
text(ax1, 0.02, 0.98, txt, 'Units','normalized', ...
     'VerticalAlignment','top','HorizontalAlignment','left', ...
     'FontName','FixedWidth','FontSize',9, 'Color','w', ...
     'BackgroundColor',[0 0 0 0.55], 'Margin',6);

% Panel B: the velocity trace that located the event
ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
nF  = numel(kin.v);
i1  = max(1,   impact - opt.Pad);
i2  = min(nF,  stopF  + opt.Pad);
idx = i1:i2;
t_ms = kin.t_s(idx)*1000;

plot(ax2, t_ms, kin.v(idx), '-', 'Color',[0.15 0.15 0.15], ...
     'LineWidth',1.4, 'DisplayName','v (smoothed)');
plot(ax2, kin.t_s(impact)*1000, kin.v(impact), 'v', 'MarkerSize',10, ...
     'MarkerFaceColor',[0.20 0.80 1.00],'MarkerEdgeColor','k', ...
     'DisplayName',sprintf('impact (%d)',impact));
if stopF >= 1 && stopF <= nF
    plot(ax2, kin.t_s(stopF)*1000, kin.v(stopF), 'o', 'MarkerSize',9, ...
        'MarkerFaceColor',[1 0.45 0.1],'MarkerEdgeColor','k', ...
        'DisplayName',sprintf('stopFrame (%d)',stopF));
end
yline(ax2, 0, ':', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');

% t_stop is the KD linear extrapolation of the pre-crossing segment, not the
% zero-crossing frame itself; draw the line so the two can be told apart.
if isfinite(local_get(kin,'a_stop_cm_s2',NaN)) && kin.a_stop_cm_s2 ~= 0
    tS  = kin.t_stop_s*1000;
    tt  = linspace(min(t_ms), max(tS + 2, max(t_ms)), 50);
    vv  = kin.a_stop_cm_s2 * (tt/1000 - kin.t_stop_s);
    plot(ax2, tt, vv, '--', 'Color',[0.85 0.33 0.10], 'LineWidth',1.3, ...
        'DisplayName','KD extrapolation  v = a\_stop (t - t\_stop)');
    xline(ax2, tS, '-.', 'Color',[0.85 0.33 0.10], ...
        'HandleVisibility','off');
    text(ax2, tS, 0, sprintf('  t\\_stop = %.3f ms', tS), ...
        'Color',[0.85 0.33 0.10], 'VerticalAlignment','bottom', 'FontSize',8);
end
ylim(ax2, local_pad_ylim(kin.v(idx)));
xlabel(ax2,'t - t_{impact}  (ms)'); ylabel(ax2,'v  (cm/s)');
title(ax2, sprintf('velocity: impact %d, stop %d  (pad %d)', impact, stopF, opt.Pad));
legend(ax2,'Location','best','Box','off');

% ── 5) Optional PNG ──────────────────────────────────────────────────────
outPath = '';
if opt.Save
    outDir = fullfile(opt.Root,'03_RESULTS','_batch_logs','impact_checks');
    if ~isfolder(outDir), mkdir(outDir); end
    outPath = fullfile(outDir, [trialTag '_impact_check.png']);
    exportgraphics(fig, outPath, 'Resolution', 150);
    fprintf('  wrote %s\n', outPath);
end

if ~opt.Show
    close(fig);
    fig = gobjects(1);      % nothing to hand back; the PNG is the output
end

fprintf('\n');
R = struct('trialTag',trialTag, 'impact_index',impact, 'stopFrame',stopF, ...
           'rawFrame',rawFrame, 'refMarkerID',refID, ...
           'rodBedDist_px_at_impact',rbAtImpact, ...
           'impactDistPx',calib.impactDistPx, ...
           'offset_px',rbAtImpact - calib.impactDistPx, ...
           'v0_cm_s',kin.v0_cm_s, 'd_final_cm',kin.d_final_cm, ...
           'calibSource',calibSrc, 'frameSource',frameSrc, ...
           'tracksPath',tracksPath, 'figurePath',outPath, 'figure',fig);
end

% ─────────────────────────────────────────────────────────────────────────
function [img, src] = local_get_frame(tracksPath, meta, trialTag, rawFrame, opt)
%LOCAL_GET_FRAME  Exported PNG if it survives, else straight from the raw avi.
img = []; src = '';

% 1) The exported PNG. framesDir mirrors the results path under 01_FRAMES, and
%    may sit under a local scratch root (JERBOA_FRAMES_ROOT), so try both.
trackDir  = fileparts(tracksPath);              % <resultsDir>/tracks
resultDir = fileparts(trackDir);                % <resultsDir>
name      = sprintf('frame_%05d.png', rawFrame);
cands     = {};
if contains(resultDir, '03_RESULTS')
    cands{end+1} = fullfile(strrep(resultDir,'03_RESULTS','01_FRAMES'), name);
    fr = getenv('JERBOA_FRAMES_ROOT');
    if ~isempty(fr)
        rel = extractAfter(string(resultDir), "03_RESULTS");
        cands{end+1} = char(fullfile(fr, '01_FRAMES', strip(rel, filesep), name));
    end
end
for k = 1:numel(cands)
    if isfile(cands{k})
        img = imread(cands{k});
        src = sprintf('exported PNG (%s)', cands{k});
        return
    end
end

% 2) Fall back to the raw video. Exact stem match including the material and
%    container folders: a substring test on the height label is not enough,
%    since '25mm' is contained in '125mm' and '325mm'.
V = dir(fullfile(opt.RawRoot,'**','*.avi'));
V = V(~[V.isdir]);
if isempty(V)
    error('diag_impact_frame:noFrameNoVideo', ...
        ['Frame %d not found under 01_FRAMES and no .avi under %s. ' ...
         'Re-export frames or point RawRoot at the capture tree.'], ...
        rawFrame, opt.RawRoot);
end
paths = unique(string(fullfile({V.folder}', {V.name}')));
[~, stems] = arrayfun(@(p) fileparts(p), paths, 'UniformOutput', false);
stems = string(stems);
stem  = sprintf('%s_T%02d', meta.heightLabel, meta.trialNum);
matDir  = string(filesep) + string(meta.material)  + string(filesep);
contDir = string(filesep) + string(meta.container) + string(filesep);
hit = strcmpi(stems, stem) & contains(paths, matDir,  'IgnoreCase', true) ...
                           & contains(paths, contDir, 'IgnoreCase', true);
if ~any(hit)
    error('diag_impact_frame:noMatch', ...
        ['Frame %d not found under 01_FRAMES, and no raw video with stem "%s" ' ...
         'under a %s%s path in %s.'], rawFrame, stem, matDir, contDir, opt.RawRoot);
end
cand = paths(hit);
orig = cand(~contains(cand, "transcoded", 'IgnoreCase', true));
if ~isempty(orig), cand = orig; end
videoPath = char(cand(1));

v = open_video(videoPath);
v.CurrentTime = (rawFrame - 1) / v.FrameRate;
if ~hasFrame(v)
    error('diag_impact_frame:frameBeyondEnd', ...
        'Raw frame %d is past the end of %s.', rawFrame, videoPath);
end
img = readFrame(v);
src = sprintf('raw video (%s)', videoPath);
fprintf('  (frames purged -- read frame straight from the .avi)\n');
end

function y = local_pad_ylim(v)
v = v(isfinite(v));
if isempty(v), y = [-1 1]; return; end
lo = min(v); hi = max(v); r = hi - lo;
if r == 0, r = max(abs(hi),1); end
y = [lo - 0.1*r, hi + 0.1*r];
end

function v = local_get(s, f, dflt)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end

function r = local_default_raw_root()
%LOCAL_DEFAULT_RAW_ROOT  Where the raw capture tree lives.
%   Set JERBOA_RAW_ROOT to point the toolchain at the current campaign without
%   editing code; the fallback is the campaign-1 tree. The lookup that uses
%   this searches recursively, so it need only be an ANCESTOR of the clips.
    r = getenv('JERBOA_RAW_ROOT');
    if isempty(r), r = 'D:\ME_GRANULAB\Test Batches'; end
end
