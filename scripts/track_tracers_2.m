%% SECTION 1: SETUP + ROD AXIS + MANUAL BED (2 clicks)
clear; close all; clc;

%% ---------------- USER IDENTIFIERS ----------------
projectRoot = ['/Users/muhannadalsharji/Library/CloudStorage/GoogleDrive-msharji@umich.edu/' ...
                'My Drive/Research /Videos/SP26/JerboaImpact_VideoExports'];
material    = 'GB';        % 'GB' or 'CHIN'
batchName   = 'Batch1';
heightLabel = 'H06';
trialNum    = 4;           % numeric
clipType    = 'long';      % 'long' or 'short'
trialID = sprintf('T%02d', trialNum);
version = lower(clipType);

baseName = sprintf('%s_%s_%s', heightLabel, trialID, version);

%% ---------------- INPUT VIDEO PATH ----------------
switch version
    case 'long'
        inDir = fullfile(projectRoot,'03_EXPORTS_LONG',heightLabel);
    case 'short'
        inDir = fullfile(projectRoot,'04_EXPORTS_SHORT',heightLabel);
    otherwise
        error('clipType must be long or short');
end

videoFile = sprintf('%s_%s_%s.mov', heightLabel, trialID, version);
videoPath = fullfile(inDir, videoFile);

if ~exist(videoPath,'file')
    error('Video not found: %s', videoPath);
end

% ---------------- RESULTS ROOT ----------------
outRoot = fullfile(projectRoot,'07_RESULTS',material,batchName);

subFolders = {'detections','tracks','kinematics','figures','qa','logs'};
for k = 1:numel(subFolders)
    f = fullfile(outRoot,subFolders{k});
    if ~exist(f,'dir'), mkdir(f); end
end

% calibration (edit if needed) (NOTE: I double checked this and it seems to
% be consistent!)
mm_per_px = 0.084682;

% ---------------- LOAD DETECTIONS ----------------
detPath = fullfile(outRoot,'detections',[baseName '_detections.mat']);
if ~exist(detPath,'file')
    error('Detections MAT not found: %s', detPath);
end

S = load(detPath);
det = S.det;

% ---------------- REFERENCE FRAME ----------------
refFrame0 = 0;
iRef = refFrame0 + 1;

v = VideoReader(videoPath);
ref_img = read(v,iRef);

figure(1); clf; imshow(ref_img); hold on;
title(sprintf('%s reference frame0=%d', baseName, refFrame0));

% ---------------- CLICK ROD AXIS ----------------
disp('Click rod TOP, then rod TIP');
[xr, yr] = ginput(2);

rod_top = [xr(1) yr(1)];
rod_tip = [xr(2) yr(2)];

plot(rod_top(1),rod_top(2),'go','LineWidth',2);
plot(rod_tip(1),rod_tip(2),'ro','LineWidth',2);
plot([rod_top(1) rod_tip(1)], [rod_top(2) rod_tip(2)], 'b-','LineWidth',2);

% rod axis unit vector
u_ref = rod_top - rod_tip;
u_ref = u_ref / norm(u_ref);

% ---------------- CLICK BED LINE (VERTICAL MODEL) ----------------
disp('Click TWO points along the bed surface (TOP then BOTTOM)');  % pick along the vertical boundary
[xb, yb] = ginput(2);

plot(xb, yb, 'y-', 'LineWidth', 2);

% Bed modeled as: x_bed(y) = m_bed_y * y + b_bed_x
m_bed_y = (xb(2) - xb(1)) / (yb(2) - yb(1) + eps);
b_bed_x = xb(1) - m_bed_y * yb(1);

fprintf('Bed: x = (%.6f)*y + (%.3f)\n', m_bed_y, b_bed_x);

% Sanity: x_bed across full height should sit inside [0..W]

H = size(ref_img,1);  W = size(ref_img,2);
y_test = [1 H];
x_test = m_bed_y*y_test + b_bed_x;
fprintf('x_bed at y=1..H: %.2f to %.2f px (should be within 0..W=%d)\n', x_test(1), x_test(2), W);

% overlay sanity (VERTICAL bed: x = m*y + b)

y_line = [1 size(ref_img,1)];          % y from top to bottom of image
x_line = m_bed_y * y_line + b_bed_x;   % x along the bed line
plot(x_line, y_line, 'y-', 'LineWidth', 2);

fprintf('m_bed_y=%.6f, b_bed_x=%.3f\n', m_bed_y, b_bed_x);
fprintf('x_bed at y=1..H: %.2f to %.2f px (should be within 0..W=%d)\n', ...
    x_line(1), x_line(2), size(ref_img,2));

disp('Section 1 is done.')
%% SECTION 2: MARKER IDS FROM REFERENCE FRAME

% --- Extract detections at reference frame ---
Cref = det.detect.centersCell{iRef};

if isempty(Cref)
    error(['No circles detected in reference frame0 = %d.\n' ...
           'Change refFrame0 and rerun Section 1.'], refFrame0);
end

% Separate x and y
x0 = Cref(:,1);
y0c = Cref(:,2);

% --- Compute rod axis unit vector (from Section 1 clicks) ---
u_ref = rod_top - rod_tip;
u_ref = u_ref / norm(u_ref);

% --- Project detections along rod axis (tip → top ordering) ---
proj = ([x0 y0c] - rod_tip) * u_ref.';
[~, idx] = sort(proj);     % smallest projection = closest to tip

x0 = x0(idx);
y0c = y0c(idx);

nMarkers = numel(x0);
markerID = (1:nMarkers)';   % ID 1 = closest to rod tip

% --- Store marker reference data ---
markers = struct();
markers.ID        = markerID;
markers.x_px0     = x0;
markers.y_px0     = y0c;
markers.x_mm0     = x0 * mm_per_px;
markers.y_mm0     = y0c * mm_per_px;
markers.colors    = lines(nMarkers);

% Store reference geometry for reproducibility
markers.u_ref     = u_ref;
markers.rod_tip0  = rod_tip;
markers.rod_top0  = rod_top;
markers.refFrame0 = refFrame0;

% --- Visual sanity check ---
figure(1); hold on;

for k = 1:nMarkers
    c = markers.colors(k,:);
    plot(x0(k), y0c(k), 'o', ...
        'MarkerSize', 8, ...
        'LineWidth', 2, ...
        'Color', c, ...
        'MarkerFaceColor', 'none');
    
    text(x0(k)+5, y0c(k), sprintf('%d', markerID(k)), ...
        'Color', c, ...
        'FontSize', 9, ...
        'FontWeight','bold');
end

disp(['Section 2 complete: assigned ', num2str(nMarkers), ' marker IDs.']);

%% SECTION 3: TRACK MARKERS (rigid 1D translation, robust)

v = VideoReader(videoPath);

nFrames = floor(v.Duration * v.FrameRate);
fps_used = v.FrameRate;

nMarkers = numel(markers.ID);

% --- Preallocate ---
x_px = nan(nMarkers, nFrames);
y_px = nan(nMarkers, nFrames);
s_px = nan(nMarkers, nFrames);

rod_tip_px = nan(nFrames,2);
z_tip_mm   = nan(nFrames,1);

% --- Reference geometry ---
s_ref = ([markers.x_px0 markers.y_px0] - markers.rod_tip0) * markers.u_ref.';
rod_tip_px_ref = markers.rod_tip0;

assignTol_loose = 120;   % only for shift estimation
assignTol_tight = 30;    % final assignment gate

delta_s_last = 0;        % predicted translation (px along rod)

for fi = 1:nFrames

    C = det.detect.centersCell{fi};
    if isempty(C)
        continue
    end

    % Project detections
    s_det = (C - rod_tip_px_ref) * markers.u_ref.';

    % ---------- PASS 1: estimate translation ----------
    usedDet = false(size(s_det));
    s_temp  = nan(nMarkers,1);

    for k = 1:nMarkers
        ds = abs(s_det - (s_ref(k) + delta_s_last));
        ds(usedDet) = inf;

        [dmin,j] = min(ds);

        if dmin < assignTol_loose
            s_temp(k) = s_det(j);
            usedDet(j) = true;
        end
    end

    valid = isfinite(s_temp);

    if any(valid)
        delta_s = median(s_temp(valid) - s_ref(valid));
        delta_s_last = delta_s;
    else
        delta_s = delta_s_last;
    end

    % ---------- PASS 2: final assignment ----------
    usedDet = false(size(s_det));

    for k = 1:nMarkers
        target = s_ref(k) + delta_s;

        ds = abs(s_det - target);
        ds(usedDet) = inf;

        [dmin,j] = min(ds);

        if dmin < assignTol_tight
            x_px(k,fi) = C(j,1);
            y_px(k,fi) = C(j,2);
            s_px(k,fi) = s_det(j);
            usedDet(j) = true;
        end
    end

    % ---------- Rod tip position ----------
    rod_tip_px(fi,:) = rod_tip_px_ref + delta_s * markers.u_ref;

    % ---------- Depth relative to bed (VERTICAL BED) ----------
    x_now = rod_tip_px(fi,1);
    y_now = rod_tip_px(fi,2);
    
    xbed = m_bed_y * y_now + b_bed_x;   % bed x at this y
    dx   = x_now - xbed;                % + means to the RIGHT of bed
    
    % Choose sign so penetration into bed is NEGATIVE (UFL-style)
    z_tip_mm(fi) = +dx * mm_per_px;

end

disp('Section 3 complete: rigid tip trajectory and z(t) computed.');

%% SECTION 4–6: t(t0=impact) + v(t) + a(t) + t_stop (UFL-style)

% SECTION 4: Time + units
impactFrame0 = 6;                       % PLEASE MAKE SURE TO EDIT PER TRIAL!!
g_cm_s2 = 980;
fps_true= 2715;                          % true camera acquisition fps 

frame0 = (0:nFrames-1)';                 % 0-based
t_s    = (frame0 - impactFrame0) / fps_true;   % aligned so t=0 at impact
impact_index = impactFrame0 + 1;         % MATLAB index

z_tip_cm = z_tip_mm(:) / 10;             % UFL-style reporting
z_valid_mask = isfinite(z_tip_cm);

%% SECTION 5: v(t) via local linear fit of z(t)
% window: 21 frames

Wz = 21;
halfWz = floor(Wz/2);

v_tip_cm_s = nan(nFrames,1);

for i = 1:nFrames
    i1 = max(1, i-halfWz);
    i2 = min(nFrames, i+halfWz);

    tt = t_s(i1:i2);
    zz = z_tip_cm(i1:i2);

    ok = isfinite(tt) & isfinite(zz);
    if sum(ok) < 3
        continue
    end

    p = polyfit(tt(ok), zz(ok), 1);   % z ≈ p1*t + p0
    v_tip_cm_s(i) = p(1);             % slope = velocity
end

% v0 at impact: local-fit velocity at the impact index
v0_cm_s = v_tip_cm_s(impact_index);

%% SECTION 6: a(t) via local linear fit of v(t)
% window: 31 frames

Wv = 31;
halfWv = floor(Wv/2);

a_tip_cm_s2 = nan(nFrames,1);

for i = 1:nFrames
    i1 = max(1, i-halfWv);
    i2 = min(nFrames, i+halfWv);

    tt = t_s(i1:i2);
    vv = v_tip_cm_s(i1:i2);

    ok = isfinite(tt) & isfinite(vv);
    if sum(ok) < 3
        continue
    end

    p = polyfit(tt(ok), vv(ok), 1);    % v ≈ p1*t + p0
    a_tip_cm_s2(i) = p(1);             % slope = acceleration
end

a_plus_g_cm_s2 = a_tip_cm_s2 + g_cm_s2;

% Stop time (official): extrapolation on velocity band
% band: 1–10 cm/s, after impact, penetrating only (v<0)

stop_band = [1 10];             % cm/s
minFitFrames = 10;

post = (t_s >= 0);
vfit = v_tip_cm_s;

use = post & isfinite(vfit) & (vfit < 0) & ...
      (abs(vfit) >= stop_band(1)) & (abs(vfit) <= stop_band(2));

t_stop_s = NaN;
stop_fit_coeff = [NaN NaN];     % [p1 p0]
stop_nFitFrames = sum(use);

if stop_nFitFrames >= minFitFrames
    p = polyfit(t_s(use), vfit(use), 1);     % v ≈ p1*t + p0
    p1 = p(1); p0 = p(2);

    % Reject nonsense fits
    if isfinite(p1) && isfinite(p0) && (p1 < 0)
        t_candidate = -p0 / p1;

        % Must land within the time span of the fit band
        tmin = min(t_s(use));
        tmax = max(t_s(use));

        if (t_candidate >= tmin) && (t_candidate <= tmax)
            t_stop_s = t_candidate;
            stop_fit_coeff = [p1 p0];
        end
    end
end

% ----------------------------
% Stop time sanity check (threshold)
% |v_fit| < 1 cm/s for 20 consecutive frames after impact
% ----------------------------
thr_cm_s = 1;
persistN = 20;

t_stop_thr_s = NaN;
stop_discrepancy_flag = false;

idxPost = find(t_s >= 0);
if ~isempty(idxPost)
    vvPost = abs(vfit(idxPost));
    okPost = isfinite(vvPost);

    below = false(size(vvPost));
    below(okPost) = (vvPost(okPost) < thr_cm_s);

    % Find first run of persistN trues
    runLen = 0;
    for j = 1:numel(below)
        if below(j)
            runLen = runLen + 1;
            if runLen >= persistN
                j0 = j - persistN + 1;
                t_stop_thr_s = t_s(idxPost(j0));
                break
            end
        else
            runLen = 0;
        end
    end
end

% Flag if disagreement > 20 ms (0.02 s)
if isfinite(t_stop_s) && isfinite(t_stop_thr_s)
    stop_discrepancy_flag = (abs(t_stop_s - t_stop_thr_s) > 0.02);
end

% Optional convenience outputs
t_stop_index = NaN;
z_at_stop_cm = NaN;

if isfinite(t_stop_s)
    [~, t_stop_index] = min(abs(t_s - t_stop_s));
    z_at_stop_cm = z_tip_cm(t_stop_index);
end

% ----------------------------
% Pack for saving to _kin.mat later
% ----------------------------
kin = struct();

kin.meta.projectRoot  = projectRoot;
kin.meta.material     = material;
kin.meta.batchName    = batchName;
kin.meta.heightLabel  = heightLabel;
kin.meta.trialID      = trialID;
kin.meta.version      = version;
kin.meta.videoPath    = videoPath;
kin.meta.detPath      = detPath;

kin.calib.mm_per_px = mm_per_px;

% bed model 
kin.bed.model = 'x = m*y + b';
kin.bed.m_bed_y = m_bed_y;
kin.bed.b_bed_x = b_bed_x;

kin.params.fps_used      = fps_used;
kin.params.impactFrame0  = impactFrame0;
kin.params.g_cm_s2       = g_cm_s2;

kin.time.frame0       = frame0;
kin.time.t_s          = t_s;
kin.time.nFrames      = nFrames;
kin.time.impact_index = impact_index;

kin.pos.rod_tip_px    = rod_tip_px;
kin.pos.z_tip_cm      = z_tip_cm;
kin.pos.z_valid_mask  = z_valid_mask;

kin.vel.v_tip_cm_s    = v_tip_cm_s;
kin.vel.v0_cm_s       = v0_cm_s;
kin.vel.v_fit_window  = Wz;
kin.vel.v_method      = 'local_linear_fit_z';

kin.acc.a_tip_cm_s2       = a_tip_cm_s2;
kin.acc.a_plus_g_cm_s2    = a_plus_g_cm_s2;
kin.acc.a_fit_window      = Wv;
kin.acc.a_method          = 'local_linear_fit_v';

kin.stop.t_stop_s            = t_stop_s;
kin.stop.stop_fit_coeff      = stop_fit_coeff;   % [p1 p0]
kin.stop.stop_fit_band_cm_s  = stop_band;
kin.stop.stop_nFitFrames     = stop_nFitFrames;

kin.stop.t_stop_thr_s           = t_stop_thr_s;
kin.stop.stop_discrepancy_flag  = stop_discrepancy_flag;

kin.stop.t_stop_index = t_stop_index;
kin.stop.z_at_stop_cm = z_at_stop_cm;

kin.track.x_px = x_px;
kin.track.y_px = y_px;
kin.track.s_px = s_px;          % why not
kin.track.nMarkers = nMarkers;
kin.track.markers = markers;    % contains IDs, reference geometry, u_ref, etc.

% --- save ---
kinDir = fullfile(outRoot,'kinematics',heightLabel);
if ~exist(kinDir,'dir'), mkdir(kinDir); end

kinPath = fullfile(kinDir, [baseName '_kin.mat']);
save(kinPath,'kin','-v7.3');
fprintf('Saved: %s\n', kinPath);

disp('Sections 4–6 complete: t, z, v, a, and stop-time computed.');

%% --- SANITY PLOTS: z, v, a (subplots) ---
% Assumes you already have:
% t_s, z_tip_cm, v_tip_cm_s, a_tip_cm_s2, impactFrame0, fps_true, nFrames

impact_index = impactFrame0 + 1;

% ----- quick prints -----
fprintf('\n--- SANITY (z,v,a) ---\n');
fprintf('t_s(impact_index) = %.12f s (should be 0)\n', t_s(impact_index));
fprintf('dt mean/std       = %.12g / %.12g s (should be 1/fps_true)\n', mean(diff(t_s)), std(diff(t_s)));
fprintf('1/fps_true         = %.12g s\n', 1/fps_true);

fprintf('z_tip_cm min/max   = %.3f / %.3f cm\n', min(z_tip_cm,[],'omitnan'), max(z_tip_cm,[],'omitnan'));
fprintf('v_tip_cm_s min/max = %.3f / %.3f cm/s\n', min(v_tip_cm_s,[],'omitnan'), max(v_tip_cm_s,[],'omitnan'));
fprintf('a_tip_cm_s2 min/max= %.1f / %.1f cm/s^2\n', min(a_tip_cm_s2,[],'omitnan'), max(a_tip_cm_s2,[],'omitnan'));

% ----- sign sanity right after impact (should match your convention) -----
k1 = min(nFrames, impact_index + 2);
k2 = min(nFrames, impact_index + 30);

dz = z_tip_cm(k2) - z_tip_cm(k1);                    % should be < 0 if penetration is negative
vmed = median(v_tip_cm_s(k1:k2), 'omitnan');         % should be < 0
dv = v_tip_cm_s(k2) - v_tip_cm_s(k1);                % should be > 0 as it slows toward 0

fprintf('post-impact dz = %.5f cm (want < 0)\n', dz);
fprintf('post-impact median v = %.5f cm/s (want < 0)\n', vmed);
fprintf('post-impact dv = %.5f cm/s (want > 0)\n', dv);
fprintf('---------------------\n\n');

% ----- re-zero z at impact (so z_rel(t=0)=0) -----
z_rel_cm = z_tip_cm - z_tip_cm(impact_index);

% ----- full timeline subplots -----
figure('Color','w');
subplot(3,1,1);
plot(t_s, z_rel_cm, '.-'); grid off; xline(0,'k--');
ylabel('z - z(0) (cm)');
title('Full timeline (re-zeroed: penetration should go negative)');

subplot(3,1,2);
plot(t_s, v_tip_cm_s, '.-'); grid off; xline(0,'k--');
ylabel('v (cm/s)');

subplot(3,1,3);
plot(t_s, a_tip_cm_s2, 'o'); grid off; xline(0,'k--');
xlabel('t (s)'); ylabel('a (cm/s^2)');

% ----- zoom around impact -----
W = 0.05;  % seconds on each side
mask = (t_s >= -W) & (t_s <= W);

figure('Color','w');
subplot(3,1,1);
plot(t_s(mask), z_rel_cm(mask), '.-'); grid off; xline(0,'k--');
ylabel('z - z(0) (cm)');
title(sprintf('Zoom around impact (\\pm %.3f s)', W));

subplot(3,1,2);
plot(t_s(mask), v_tip_cm_s(mask), '.-'); grid off; xline(0,'k--');
ylabel('v (cm/s)');

subplot(3,1,3);
plot(t_s(mask), a_tip_cm_s2(mask), '.-'); grid off; xline(0,'k--');
xlabel('t (s)'); ylabel('a (cm/s^2)');

%% ----- Sanity Check Prompts -----%%

% A) Bed geometry sanity; Expect: bed overlay matches the physical wall; 
% x_bed stays within 0..W; |Δx| small (nearly vertical)

H = size(ref_img,1);  W = size(ref_img,2);

y_line = [1 H];
x_line = m_bed_y*y_line + b_bed_x;

fprintf('\n--- BED SANITY ---\n');
fprintf('Frame W=%d, H=%d\n', W, H);
fprintf('Bed model: x = m*y + b  (m=%.6g, b=%.3f)\n', m_bed_y, b_bed_x);
fprintf('x_bed at y=1..H: %.2f -> %.2f px (should be within 0..%d)\n', x_line(1), x_line(2), W);
fprintf('|Δx| across height: %.3f px (vertical wall => small)\n', abs(x_line(2)-x_line(1)));

figure('Color','w'); imshow(ref_img); hold on;
plot(x_line, y_line, 'y-', 'LineWidth', 2);
title('Bed overlay (yellow) on reference frame');

%% B) Penetration sign sanity (MUST pass); Expect: post-impact dz<0 and post-impact median v<0 

impact_index = impactFrame0 + 1;
k1 = impact_index + 2;
k2 = min(nFrames, impact_index + 30);

dz = z_tip_cm(k2) - z_tip_cm(k1);
vmed = median(v_tip_cm_s(k1:k2), 'omitnan');

fprintf('\n--- SIGN SANITY ---\n');
fprintf('Post-impact dz = %.6f cm (want < 0)\n', dz);
fprintf('Post-impact median v = %.6f cm/s (want < 0)\n', vmed);

%% C) Pixel→mm sanity using rod length: Expect: measured pixel distance 
% between your two clicks ≈ rod_length_mm / mm_per_px. (should be less than
% 2 or 5% 

rod_length_mm = 28.15; % your known rod length

rod_len_px = norm(rod_top - rod_tip);
mm_per_px_meas = rod_length_mm / rod_len_px;

fprintf('\n--- mm/px SANITY (rod length) ---\n');
fprintf('rod_len_px = %.3f px\n', rod_len_px);
fprintf('mm_per_px (current) = %.6f\n', mm_per_px);
fprintf('mm_per_px (measured from clicks) = %.6f\n', mm_per_px_meas);
fprintf('percent diff = %.2f %%\n', 100*(mm_per_px_meas - mm_per_px)/mm_per_px);

%% D) Tracking completeness sanity: Expect: most frames have many assigned markers; 
% if median is like 4/8 consistently, tracking is shaky

nAssigned = sum(isfinite(x_px) & isfinite(y_px), 1);
fprintf('\n--- TRACKING SANITY ---\n');
fprintf('assigned per frame min/median/max = %d / %d / %d (out of %d)\n', ...
    min(nAssigned), round(median(nAssigned,'omitnan')), max(nAssigned), size(x_px,1));

figure('Color','w'); plot(nAssigned,'.-'); grid on;
xlabel('frame'); ylabel('markers assigned'); title('Markers assigned per frame');

%% E) Dynamics plausibility sanity

fprintf('\n--- DYNAMICS RANGE ---\n');
fprintf('z_tip_cm min/max   = %.3f / %.3f cm\n', min(z_tip_cm,[],'omitnan'), max(z_tip_cm,[],'omitnan'));
fprintf('v_tip_cm_s min/max = %.3f / %.3f cm/s\n', min(v_tip_cm_s,[],'omitnan'), max(v_tip_cm_s,[],'omitnan'));
fprintf('a_tip_cm_s2 min/max= %.1f / %.1f cm/s^2\n', min(a_tip_cm_s2,[],'omitnan'), max(a_tip_cm_s2,[],'omitnan'));

% Post-impact "settled" check: last 200 frames should be near constant z and v ~ 0
tail = max(1,nFrames-200):nFrames;
fprintf('tail z std = %.5f cm (small is good)\n', std(z_tip_cm(tail),'omitnan'));
fprintf('tail v median = %.5f cm/s (near 0 is good)\n', median(v_tip_cm_s(tail),'omitnan'));
