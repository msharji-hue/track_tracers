%% ========================================================================
%  PER-MODEL SCALE / TRACKING DIAGNOSTIC  (GB/full, Tight | Default | Wide)
%
%  The rod's inter-marker spacing is a FIXED PHYSICAL LENGTH. Measuring it in
%  pixels per model therefore calibrates the scale directly:
%
%      mmPerPx(model) = mmPerPx(ref) * spacing_px(ref) / spacing_px(model)
%
%  It also tests the swap hypothesis: trackTolerancePx must stay well BELOW one
%  spacing, or a stale reference can grab a neighbouring marker and inject a
%  one-spacing step into z (which reads as a fake impact + velocity spike).
%
%  Read-only. Writes nothing.
%% ========================================================================

root   = 'D:\ME_GRANULAB\JerboaImpact';
gbRoot = fullfile(root,'03_RESULTS','GB','Batch 5');
base   = get_calibration();

D = dir(fullfile(gbRoot,'**','*_tracks.mat'));
n = numel(D);
mdl=strings(n,1); tag=mdl; spac=nan(n,1); span=nan(n,1); nmk=nan(n,1);
fpsv=nan(n,1); hmm=nan(n,1); maxJump=nan(n,1);

for k = 1:n
    tp = fullfile(D(k).folder, D(k).name);
    try, S = load(tp,'meta','tracks'); catch, continue; end
    m = S.meta;
    if ~strcmpi(getf(m,'container',''),'full'), continue; end

    tag(k) = string(getf(m,'trialTag',''));
    mdl(k) = model_of(tag(k), string(D(k).folder));
    hmm(k) = getf(m,'dropHeight_mm',NaN);
    fpsv(k)= getf(m,'fps_true',NaN);

    X = S.tracks.trackedX; Y = S.tracks.trackedY;
    bn = sqrt(base.lineA^2 + base.lineB^2);
    d  = (base.lineA.*X + base.lineB.*Y + base.lineC)./bn;

    ff = find(all(isfinite(X),1),1,'first');
    if isempty(ff), continue; end
    col = sort(d(:,ff));
    nmk(k)  = numel(col);
    spac(k) = median(abs(diff(col)));      % px between adjacent markers
    span(k) = max(col) - min(col);

    % largest single-frame jump of the reference marker (swap detector)
    colf = d(:,ff); colf(~isfinite(colf)) = +Inf;
    [~,ref] = min(colf);
    r = d(ref,:); r = r(isfinite(r));
    if numel(r) > 2, maxJump(k) = max(abs(diff(r))); end
end

ok = tag ~= "";
Tm = table(mdl(ok),tag(ok),hmm(ok),fpsv(ok),nmk(ok),spac(ok),span(ok),maxJump(ok), ...
    'VariableNames',{'model','trialTag','dropHeight_mm','fps','nMarkers', ...
                     'spacing_px','span_px','maxFrameJump_px'});

%% ---- 1. spacing per model ---------------------------------------------
fprintf('\n========== MARKER SPACING BY MODEL ==========\n');
G = groupsummary(Tm,'model',{'median','mean','std'},{'spacing_px','span_px'});
disp(G)

ref = "Default";
sRef = median(Tm.spacing_px(Tm.model==ref),'omitnan');
fprintf('reference: %s spacing = %.2f px  at mmPerPx = %.4f\n', ref, sRef, base.mmPerPx);
fprintf('  -> physical spacing = %.3f mm\n\n', sRef*base.mmPerPx);

fprintf('IMPLIED mmPerPx per model (rod spacing is a fixed length):\n');
for mm_ = unique(Tm.model)'
    s = median(Tm.spacing_px(Tm.model==mm_),'omitnan');
    implied = base.mmPerPx * sRef / s;
    fprintf('  %-8s spacing %6.2f px  ->  mmPerPx = %.4f   (%.1f%% vs %.4f)\n', ...
        mm_, s, implied, 100*(implied/base.mmPerPx - 1), base.mmPerPx);
end

%% ---- 2. swap risk ------------------------------------------------------
fprintf('\n========== SWAP RISK ==========\n');
fprintf('trackTolerancePx = %g\n', base.trackTolerancePx);
fprintf('A swap becomes possible when tolerance approaches one spacing.\n\n');
for mm_ = unique(Tm.model)'
    s = median(Tm.spacing_px(Tm.model==mm_),'omitnan');
    j = Tm.maxFrameJump_px(Tm.model==mm_);
    nSwap = sum(j > 0.6*s);                 % jump > 60% of a spacing
    fprintf('  %-8s spacing %6.2f px | tol/spacing = %.2f | trials with a jump >0.6*spacing: %d of %d\n', ...
        mm_, s, base.trackTolerancePx/s, nSwap, sum(Tm.model==mm_));
end

%% ---- 3. distribution of single-frame jumps ----------------------------
figure('Color','w','Name','single-frame jumps by model');
models = unique(Tm.model); mc = lines(numel(models));
hold on; grid on;
for i = 1:numel(models)
    j = Tm.maxFrameJump_px(Tm.model==models(i));
    j = j(isfinite(j));
    histogram(j, 40, 'FaceColor', mc(i,:), 'FaceAlpha',.5, 'DisplayName',char(models(i)));
end
s = median(Tm.spacing_px,'omitnan');
xline(s,'k--','one spacing','LineWidth',1.5,'HandleVisibility','off');
xline(base.trackTolerancePx,'r--','trackTolerancePx','HandleVisibility','off');
xlabel('largest single-frame jump of the reference marker (px)');
ylabel('trials'); legend('Location','northeast');
title('A cluster at one spacing indicates marker swaps, not motion');

%% ---- 4. spacing vs fps (did magnification change with frame rate?) ----
figure('Color','w','Name','spacing vs fps'); hold on; grid on;
for i = 1:numel(models)
    s2 = Tm(Tm.model==models(i),:);
    scatter(s2.fps, s2.spacing_px, 30, mc(i,:), 'filled', ...
            'MarkerFaceAlpha',.6,'DisplayName',char(models(i)));
end
xlabel('fps'); ylabel('inter-marker spacing (px)');
title('If spacing tracks fps, the camera was cropped/rezoomed between runs');
legend('Location','best');

fprintf('\nRead-only — nothing written.\n');

%% ---- helpers -----------------------------------------------------------
function m = model_of(tag, folder)
    t = lower(string(tag)); f = lower(string(folder));
    if endsWith(t,"_tight") || contains(f,"tight"), m = "Tight";
    elseif endsWith(t,"_wide") || contains(f,"wide"), m = "Wide";
    else, m = "Default";
    end
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
