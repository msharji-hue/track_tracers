function results = track_tracers(params)
%% track_tracers(params)
% Run via scripts/run_track_tracers.m.
% This function builds paths, loads detections, tracks markers, computes toe depth z(t),
% computes v(t), a(t), and tstop, then saves a standardized MAT output.
%
% Required helper: local_rigid2D.m must be on the MATLAB path.

if nargin < 1 || isempty(params), params = struct(); end

%% =========================
% DEFAULTS
%% =========================
def = struct();

% Paths and trial IDs
def.projectRoot = "/Users/muhannadalsharji/Library/CloudStorage/GoogleDrive-msharji@umich.edu/My Drive/Research/Videos/SP26/JerboaImpact_VideoExports";
def.material    = "GB";     % "GB" or "CHIN"
def.batchName   = "Batch1";
def.heightLabel = "H06";
def.trialNum    = 1;
def.version     = "short";  % "long" or "short"

% Physics bookkeeping
def.fps_used     = 2715;
def.impactFrame0 = 97;

% Geometry constants (locked)
def.mm_per_px           = 0.0869;
def.rod_length_mm       = 30.00;
def.marker_pitch_mm     = 4.00;
def.marker_diameter_mm  = 2.00;
def.rod_top_to_foot_mm  = 32.017;
def.rod_top_to_toe_mm   = 37.448;

% Tracking settings
def.assignGate_px      = 25;
def.minMarkersForPose  = 3;
def.useLastPoseIfBad   = true; %#ok<NASGU> (kept for future)
def.fitWinFrames       = 11;

% Click behavior
def.doClicks = true;   % if false, must provide rod/bed below

% If doClicks == false, you must provide:
% params.rod_top_px_ref = [x y];
% params.rod_tip_px_ref = [x y];
% params.m_bed          = scalar;
% params.b_bed          = scalar;

% Saving / QA
def.doSaveQA  = false;
def.doOverlay = false; %#ok<NASGU> (overlay generation not implemented here)

% Fill missing fields
fn = fieldnames(def);
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(params, f) || isempty(params.(f))
        params.(f) = def.(f);
    end
end

% Aliases
projectRoot  = char(params.projectRoot);
material     = char(params.material);
batchName    = char(params.batchName);
heightLabel  = char(params.heightLabel);
trialNum     = params.trialNum;
version      = char(params.version);

fps_used     = params.fps_used;
impactFrame0 = params.impactFrame0;

mm_per_px           = params.mm_per_px;
rod_length_mm       = params.rod_length_mm;
marker_pitch_mm     = params.marker_pitch_mm;
marker_diameter_mm  = params.marker_diameter_mm;
rod_top_to_foot_mm  = params.rod_top_to_foot_mm;
rod_top_to_toe_mm   = params.rod_top_to_toe_mm;

assignGate_px     = params.assignGate_px;
minMarkersForPose = params.minMarkersForPose;

fitWinFrames = params.fitWinFrames;
fitWinFrames = max(3, fitWinFrames);
if mod(fitWinFrames,2)==0, fitWinFrames = fitWinFrames+1; end
halfW = floor(fitWinFrames/2);

doClicks = params.doClicks;
doSaveQA = params.doSaveQA;

% Derived
toe_offset_from_rod_end_mm = rod_top_to_toe_mm - rod_length_mm;
toe_offset_px = toe_offset_from_rod_end_mm / mm_per_px;

%% =========================
% SECTION 1: PATHS + LOAD DETECTIONS + REF CLICKS (optional)
%% =========================
trialID = sprintf('T%02d', trialNum);

exportsLongDir  = fullfile(projectRoot, '03_EXPORTS_LONG',  heightLabel);
exportsShortDir = fullfile(projectRoot, '04_EXPORTS_SHORT', heightLabel);

switch lower(version)
    case 'long'
        inDir = exportsLongDir;
    case 'short'
        inDir = exportsShortDir;
    otherwise
        error('version must be ''long'' or ''short''');
end

videoFile = sprintf('%s_%s_%s.mov', heightLabel, trialID, lower(version));
videoPath = fullfile(inDir, videoFile);
if ~exist(videoPath,'file')
    error('Video not found: %s', videoPath);
end
fprintf('VideoPath: %s\n', videoPath);

resultsRoot = fullfile(projectRoot, '07_RESULTS', material, batchName);

detDir  = fullfile(resultsRoot, 'detections', heightLabel);
detPath = fullfile(detDir, sprintf('%s_%s_%s_detections.csv', heightLabel, trialID, lower(version)));
if ~exist(detPath,'file')
    error('Detections CSV not found: %s', detPath);
end
fprintf('DetPath: %s\n', detPath);

outDir = fullfile(resultsRoot, heightLabel);
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('OutDir: %s\n', outDir);

matPath     = fullfile(outDir, sprintf('%s_%s_%s_track.mat',  heightLabel, trialID, lower(version)));
qaPath      = fullfile(outDir, sprintf('%s_%s_%s_QA.png',     heightLabel, trialID, lower(version)));
overlayPath = fullfile(outDir, sprintf('%s_%s_%s_overlay.mp4',heightLabel, trialID, lower(version))); %#ok<NASGU>

det = readtable(detPath);
detDet = det(~isnan(det.x), :);  % keep only detected rows

% REF frame
v       = VideoReader(videoPath);
ref_img = read(v, 1);

if doClicks
    figure(1); clf; imshow(ref_img); hold on;
    title('REF frame (MATLAB index = 1, frame0 = 0)');

    disp('Click rod TOP, then rod TIP (rod end closest to foot, NOT toe)');
    [xr, yr] = ginput(2);
    rod_top_px_ref = [xr(1) yr(1)];
    rod_tip_px_ref = [xr(2) yr(2)];

    plot(rod_top_px_ref(1), rod_top_px_ref(2), 'go', 'MarkerSize', 8, 'LineWidth', 2);
    plot(rod_tip_px_ref(1), rod_tip_px_ref(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    plot([rod_top_px_ref(1) rod_tip_px_ref(1)], [rod_top_px_ref(2) rod_tip_px_ref(2)], 'b-', 'LineWidth', 2);

    % unit axis u_ref points TIP -> TOP
    u_ref = (rod_top_px_ref - rod_tip_px_ref);
    u_ref = u_ref / norm(u_ref);

    % toe preview
    toe_px_ref = rod_tip_px_ref - toe_offset_px * u_ref;
    plot(toe_px_ref(1), toe_px_ref(2), 'ms', 'MarkerSize', 8, 'LineWidth', 2);
    text(toe_px_ref(1)+6, toe_px_ref(2), 'toe (inferred)', 'Color','m', 'FontSize', 9);

    disp('Click two points along the bed surface line');
    [xb, yb] = ginput(2);
    plot([xb(1) xb(2)], [yb(1) yb(2)], 'y-', 'LineWidth', 2);

    m_bed = (yb(2) - yb(1)) / (xb(2) - xb(1));
    b_bed = yb(1) - m_bed * xb(1);

    W = size(ref_img,2);
    xline = [1 W];
    yline = m_bed*xline + b_bed;
    plot(xline, yline, 'y--', 'LineWidth', 1.5);

else
    if ~isfield(params,'rod_top_px_ref') || ~isfield(params,'rod_tip_px_ref') || ~isfield(params,'m_bed') || ~isfield(params,'b_bed')
        error('doClicks=false requires params.rod_top_px_ref, params.rod_tip_px_ref, params.m_bed, params.b_bed');
    end
    rod_top_px_ref = params.rod_top_px_ref;
    rod_tip_px_ref = params.rod_tip_px_ref;

    u_ref = (rod_top_px_ref - rod_tip_px_ref);
    u_ref = u_ref / norm(u_ref);

    m_bed = params.m_bed;
    b_bed = params.b_bed;
end

disp('Section 1 done.');

%% =========================
% SECTION 2: MARKER IDS FROM REF (frame0 = 0)
%% =========================
ref_rows = detDet.frame0 == 0;
x0 = detDet.x(ref_rows);
y0 = detDet.y(ref_rows);

valid = ~isnan(x0) & ~isnan(y0);
x0 = x0(valid);
y0 = y0(valid);

if isempty(x0)
    error('No circles found in REF frame (frame0 = 0).');
end

P0 = [x0 y0];
s0 = (P0 - rod_tip_px_ref) * u_ref.';   % projection along TIP->TOP axis
[~, idx] = sort(s0);                    % small s0 near tip
P0 = P0(idx,:);

nMarkers = size(P0,1);
markerID = (1:nMarkers)';

markers = struct();
markers.ID            = markerID;
markers.ref_px        = P0;                 % Nx2
markers.ref_mm        = P0 * mm_per_px;     % Nx2
markers.colors        = lines(nMarkers);
markers.u_ref         = u_ref;
markers.rod_top_px_ref= rod_top_px_ref;
markers.rod_tip_px_ref= rod_tip_px_ref;

if nMarkers ~= 8
    warning('REF detected %d markers (expected ~8). Continue, but tracking may be harder.', nMarkers);
end

if doClicks
    figure(1); hold on;
    for k = 1:nMarkers
        c = markers.colors(k,:);
        plot(P0(k,1), P0(k,2), 'o', 'MarkerSize', 8, 'LineWidth', 2, ...
            'Color', c, 'MarkerFaceColor', 'none');
        text(P0(k,1)+5, P0(k,2), sprintf('%d', markerID(k)), ...
            'Color', c, 'FontSize', 9, 'FontWeight','bold');
    end
end

disp('Section 2 done.');

%% =========================
% SECTION 3: TRACK MARKERS (fixed IDs, allow dropouts)
%% =========================
detUse = detDet;

u_ref          = markers.u_ref;       % TIP -> TOP
rod_tip_px_ref = markers.rod_tip_px_ref;

s_ref = (markers.ref_px - rod_tip_px_ref) * u_ref.';   % nMarkers x 1

frames = unique(detUse.frame0);
nF     = numel(frames);

x_px  = nan(nMarkers, nF);
y_px  = nan(nMarkers, nF);
s_px  = nan(nMarkers, nF);

nDetections    = zeros(1, nF);
nAssigned      = zeros(1, nF);
assignCostMean = nan(1, nF);
assignCostMax  = nan(1, nF);

for fi = 1:nF
    frame0_now = frames(fi);
    rows   = detUse.frame0 == frame0_now;

    xx = detUse.x(rows);
    yy = detUse.y(rows);

    valid = ~isnan(xx) & ~isnan(yy);
    xx = xx(valid);
    yy = yy(valid);

    Nd = numel(xx);
    nDetections(fi) = Nd;

    if Nd == 0
        continue
    end

    P_det = [xx yy];
    s_det = (P_det - rod_tip_px_ref) * u_ref.';     % Nd x 1

    % Predict s for each ID
    if fi == 1
        s_pred = s_ref;
    else
        s_pred = nan(nMarkers,1);
        for k = 1:nMarkers
            prevIdx = find(~isnan(s_px(k,1:fi-1)), 1, 'last');
            if isempty(prevIdx)
                s_pred(k) = s_ref(k);
            else
                s_pred(k) = s_px(k, prevIdx);
            end
        end
    end

    % Greedy assignment by |s_det - s_pred|
    [IDgrid, Dgrid] = ndgrid(1:nMarkers, 1:Nd);
    costMat = abs(s_pred(IDgrid) - s_det(Dgrid));

    [costSorted, linIdx] = sort(costMat(:), 'ascend');

    usedID  = false(nMarkers,1);
    usedDet = false(Nd,1);
    costsAccepted = [];

    for ii = 1:numel(linIdx)
        c = costSorted(ii);
        if c > assignGate_px
            break
        end

        [k, j] = ind2sub(size(costMat), linIdx(ii));
        if usedID(k) || usedDet(j)
            continue
        end

        usedID(k)  = true;
        usedDet(j) = true;

        x_px(k,fi) = P_det(j,1);
        y_px(k,fi) = P_det(j,2);
        s_px(k,fi) = s_det(j);

        costsAccepted(end+1,1) = c; %#ok<SAGROW>
    end

    nAssigned(fi) = sum(~isnan(s_px(:,fi)));

    if ~isempty(costsAccepted)
        assignCostMean(fi) = mean(costsAccepted);
        assignCostMax(fi)  = max(costsAccepted);
    end
end

x_mm = x_px * mm_per_px;
y_mm = y_px * mm_per_px;

disp('Section 3 done.');

%% =========================
% SECTION 4: POSE + TOE + z(t) + t
%% =========================
frame0 = frames(:);                           % 0-based
t      = (frame0 - impactFrame0) / fps_used;  % seconds
iImpact = find(frame0 == impactFrame0, 1, 'first'); %#ok<NASGU>

X0_all         = markers.ref_px;
rod_tip_px_ref = markers.rod_tip_px_ref;
rod_top_px_ref = markers.rod_top_px_ref;

poseValid  = false(1,nF);
R_all      = nan(2,2,nF);
tvec_all   = nan(2,nF);

rod_tip_px = nan(nF,2);
rod_top_px = nan(nF,2);
u_t        = nan(nF,2);

toe_px     = nan(nF,2);
z_toe_mm   = nan(1,nF);
z_tip_mm   = nan(1,nF);

R_last    = eye(2);
tvec_last = [0;0];

for fi = 1:nF

    vis = ~isnan(x_px(:,fi)) & ~isnan(y_px(:,fi));
    nv  = sum(vis);

    R    = R_last;
    tvec = tvec_last;

    if nv >= minMarkersForPose
        X0 = X0_all(vis,:);
        X  = [x_px(vis,fi), y_px(vis,fi)];

        [R_try, tvec_try] = local_rigid2D(X0, X);

        if all(isfinite(R_try(:))) && all(isfinite(tvec_try(:)))
            R    = R_try;
            tvec = tvec_try;

            poseValid(fi)  = true;
            R_all(:,:,fi)  = R;
            tvec_all(:,fi) = tvec;

            R_last    = R;
            tvec_last = tvec;
        end
    end

    % Translation-only endpoint update (as in your original)
    vis2 = isfinite(x_px(:,fi)) & isfinite(y_px(:,fi));
    if ~any(vis2)
        continue
    end

    X  = [x_px(vis2,fi), y_px(vis2,fi)];
    X0 = X0_all(vis2,:);

    dxy = mean(X - X0, 1);

    rod_tip_px(fi,:) = rod_tip_px_ref + dxy;
    rod_top_px(fi,:) = rod_top_px_ref + dxy;

    uvec = rod_top_px(fi,:) - rod_tip_px(fi,:);
    uvec = uvec / norm(uvec);
    u_t(fi,:) = uvec;

    toe_px(fi,:) = rod_tip_px(fi,:) - toe_offset_px * u_t(fi,:);

    % Bed line depth (your original sign convention)
    x_toe = toe_px(fi,1);  y_toe = toe_px(fi,2);
    ybed_toe = m_bed * x_toe + b_bed;

    x_tip = rod_tip_px(fi,1); y_tip = rod_tip_px(fi,2);
    ybed_tip = m_bed * x_tip + b_bed;

    z_toe_mm(fi) = -(y_toe - ybed_toe) * mm_per_px;
    z_tip_mm(fi) = -(y_tip - ybed_tip) * mm_per_px;
end

disp('Section 4 done.');

%% =========================
% SECTION 5: v(t), a(t) (UFL style) in cm/s, cm/s^2
%% =========================
z_toe_cm = z_toe_mm / 10;
t_s      = t(:).';  % row

v_toe_cm_s = nan(size(z_toe_cm));
for i = 1:(nF-1)
    if any(isnan([z_toe_cm(i), z_toe_cm(i+1), t_s(i), t_s(i+1)]))
        continue
    end
    dt = t_s(i+1) - t_s(i);
    if dt <= 0, continue, end
    v_toe_cm_s(i) = (z_toe_cm(i+1) - z_toe_cm(i)) / dt;
end
if nF >= 2
    v_toe_cm_s(end) = v_toe_cm_s(end-1);
end

a_toe_cm_s2 = nan(size(v_toe_cm_s));
v_fit_cm_s  = nan(size(v_toe_cm_s));

for i = 1:nF
    i1 = max(1, i-halfW);
    i2 = min(nF, i+halfW);

    tt = t_s(i1:i2);
    vv = v_toe_cm_s(i1:i2);

    ok2 = isfinite(tt) & isfinite(vv);
    if sum(ok2) < 3
        continue
    end

    p = polyfit(tt(ok2), vv(ok2), 1);
    a_toe_cm_s2(i) = p(1);
    v_fit_cm_s(i)  = polyval(p, t_s(i));
end

disp('Section 5 done.');

%% =========================
% SECTION 6: tstop (UFL-style) + QA plot
%% =========================
z_cm = z_toe_mm(:).' / 10;

ok = isfinite(t(:).') & isfinite(z_cm);
t_ok  = t(ok);
z_ok  = z_cm(ok); %#ok<NASGU>
f0_ok = frame0(ok); %#ok<NASGU>

v_cm_s   = nan(size(z_cm));
a_cm_s2  = nan(size(z_cm));
a_plus_g_cm_s2 = nan(size(z_cm));
tstop_s  = NaN;
tstopFrame0 = NaN;

if numel(t_ok) < 5
    warning('Not enough valid z(t) points for derivatives/tstop.');
else
    % Local linear fit on z(t) to get v(t) at valid samples
    v_ok = nan(size(t_ok));
    for ii = 1:numel(t_ok)
        iL = max(1, ii-halfW);
        iR = min(numel(t_ok), ii+halfW);
        if (iR - iL + 1) < 2, continue, end
        tt = t_ok(iL:iR);
        zz = z_cm(ok); %#ok<NASGU>  % keep consistent with original pipeline
        zz = z_ok(iL:iR);
        p = polyfit(tt, zz, 1);
        v_ok(ii) = p(1);
    end
    v_cm_s(ok) = v_ok;

    % Local linear fit on v(t) to get a(t)
    a_ok = nan(size(t_ok));
    for ii = 1:numel(t_ok)
        iL = max(1, ii-halfW);
        iR = min(numel(t_ok), ii+halfW);

        vv = v_ok(iL:iR);
        tt = t_ok(iL:iR);

        good = isfinite(vv) & isfinite(tt);
        if nnz(good) < 2, continue, end

        p = polyfit(tt(good), vv(good), 1);
        a_ok(ii) = p(1);
    end
    a_cm_s2(ok) = a_ok;

    g_cm_s2 = 981;
    a_plus_g_cm_s2 = a_cm_s2 + g_cm_s2;

    % tstop by linear extrapolation near rest
    post = isfinite(v_cm_s) & isfinite(t(:).') & (t(:).' >= 0);
    if any(post)
        vr = max(v_cm_s(post));
        if ~isfinite(vr), vr = 0; end

        sel = post & isfinite(v_cm_s) & (v_cm_s < -2*vr);
        selIdx = find(sel);

        if numel(selIdx) < 5
            selIdx = find(post & isfinite(v_cm_s));
            selIdx = selIdx(max(1, end-20):end);
        else
            selIdx = selIdx(max(1, end-20):end);
        end

        tt_fit = t(selIdx);
        vv_fit = v_cm_s(selIdx);

        tt_fit = tt_fit(:);
        vv_fit = vv_fit(:);

        good = isfinite(tt_fit) & isfinite(vv_fit);
        tt_fit = tt_fit(good);
        vv_fit = vv_fit(good);

        if numel(tt_fit) >= 2
            p = polyfit(tt_fit, vv_fit, 1);
            m = p(1); b = p(2);
            if isfinite(m) && abs(m) > 0
                tstop_s = -b/m;
                tstopFrame0 = round(tstop_s*fps_used + impactFrame0);
            end
        end
    end
end

figure(101); clf;
subplot(3,1,1);
plot(t, z_cm, 'k-'); grid on;
ylabel('z (cm)'); title('z(t), v(t), a(t) in cm units');

subplot(3,1,2);
plot(t, v_cm_s, 'r-'); grid on;
ylabel('v (cm/s)');

subplot(3,1,3);
plot(t, a_plus_g_cm_s2, 'b-'); grid on;
ylabel('a+g (cm/s^2)'); xlabel('t (s)');

disp('Section 6 done.');

%% =========================
% SECTION 7: SAVE (standardized)
%% =========================
save(matPath, ...
    'videoPath','detPath','projectRoot','resultsRoot','material','batchName','heightLabel','trialID','version', ...
    'mm_per_px','fps_used','impactFrame0','frame0','t', ...
    'm_bed','b_bed','toe_offset_from_rod_end_mm','toe_offset_px', ...
    'rod_length_mm','marker_pitch_mm','marker_diameter_mm', ...
    'rod_top_to_foot_mm','rod_top_to_toe_mm', ...
    'markers','x_px','y_px','x_mm','y_mm','s_px', ...
    'poseValid','rod_tip_px','rod_top_px','toe_px','u_t', ...
    'z_cm','v_cm_s','a_cm_s2','a_plus_g_cm_s2', ...
    'tstop_s','tstopFrame0', ...
    'nDetections','nAssigned','assignCostMean','assignCostMax');

fprintf('Saved MAT: %s\n', matPath);

if doSaveQA
    if ishandle(101)
        exportgraphics(figure(101), qaPath, 'Resolution', 300);
        fprintf('Saved QA PNG: %s\n', qaPath);
    else
        warning('Figure 101 not found. Skipping QA PNG.');
    end
end

disp('Section 7 done.');

%% =========================
% RETURN RESULTS
%% =========================
results = struct();
results.videoPath   = videoPath;
results.detPath     = detPath;
results.projectRoot = projectRoot;
results.resultsRoot = resultsRoot;
results.params      = params;

results.frame0      = frame0;
results.t           = t;

results.m_bed       = m_bed;
results.b_bed       = b_bed;

results.markers     = markers;
results.x_px        = x_px;
results.y_px        = y_px;
results.x_mm        = x_mm;
results.y_mm        = y_mm;
results.s_px        = s_px;

results.poseValid   = poseValid;
results.rod_tip_px  = rod_tip_px;
results.rod_top_px  = rod_top_px;
results.toe_px      = toe_px;
results.u_t         = u_t;

results.z_cm        = z_cm;
results.v_cm_s      = v_cm_s;
results.a_cm_s2     = a_cm_s2;
results.a_plus_g_cm_s2 = a_plus_g_cm_s2;
results.tstop_s     = tstop_s;
results.tstopFrame0 = tstopFrame0;

results.nDetections = nDetections;
results.nAssigned   = nAssigned;
results.assignCostMean = assignCostMean;
results.assignCostMax  = assignCostMax;

end