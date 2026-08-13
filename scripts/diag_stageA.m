%% ========================================================================
%  STAGE A DISSECTION — raw detections vs tracked output
%
%  The detections CSV holds every circle found per frame (before association).
%  The tracks .mat holds what association kept. Comparing them frame by frame
%  shows exactly where the motion is lost:
%
%    detections show motion + tracks show motion   -> Stage A is fine
%    detections show motion + tracks are flat      -> ASSOCIATION is the failure
%    detections are flat too                       -> DETECTION is the failure
%
%  Read-only. Writes nothing.
%% ========================================================================

root    = 'D:\ME_GRANULAB\JerboaImpact';
detRoot = fullfile(root,'02_SAVED_DETECTIONS');
gbRoot  = fullfile(root,'03_RESULTS','GB','Batch 5');
base    = get_calibration();

TAGS = ["65mm_T03_full_tight", "105mm_T08_full_tight", ...
        "65mm_T03_full_wide",  "145mm_T08_full"];        % last = Default control

for tg = TAGS
    fprintf('\n============================================================\n');
    fprintf(' %s\n', tg);
    fprintf('============================================================\n');

    % ---- raw detections -------------------------------------------------
    Dc = dir(fullfile(detRoot,'**',char(tg)+"_detections.csv"));
    if isempty(Dc)
        Dc = dir(fullfile(root,'**',char(tg)+"_detections.csv"));
    end
    if isempty(Dc)
        fprintf('  detections CSV not found — skipping\n'); continue;
    end
    Traw = readtable(fullfile(Dc(1).folder, Dc(1).name));
    fr   = Traw.frame0 + 1;                       % 1-based
    okd  = isfinite(Traw.x);

    nFdet = max(fr);
    detN  = accumarray(fr, okd, [nFdet 1]);       % detections per frame
    detXmin = accumarray(fr(okd), Traw.x(okd), [nFdet 1], @min, NaN);
    detXmax = accumarray(fr(okd), Traw.x(okd), [nFdet 1], @max, NaN);

    % ---- tracked output -------------------------------------------------
    Dt = dir(fullfile(gbRoot,'**',char(tg)+"_tracks.mat"));
    if isempty(Dt), fprintf('  tracks .mat not found — skipping\n'); continue; end
    S = load(fullfile(Dt(1).folder,Dt(1).name),'tracks','meta');
    X = S.tracks.trackedX; Y = S.tracks.trackedY;
    [nM, nFtr] = size(X);
    trkN = sum(isfinite(X),1);                    % tracked markers per frame

    fps = S.meta.fps_true;
    fvf = 1;
    if isfield(S.meta,'firstValidFrame') && isfinite(S.meta.firstValidFrame)
        fvf = S.meta.firstValidFrame;
    end

    fprintf('  detections : %d frames, %d with >=1 circle, max %d circles\n', ...
        nFdet, sum(detN>0), max(detN));
    fprintf('  tracks     : %d frames x %d markers, %d with >=1 tracked\n', ...
        nFtr, nM, sum(trkN>0));
    fprintf('  meta.firstValidFrame = %d   (tracks start here in video frames)\n', fvf);

    % ---- THE COMPARISON: motion present in each ------------------------
    dx_det = max(detXmax) - min(detXmin);
    dx_trk = NaN;
    xr = nan(nM,1);
    for j = 1:nM
        v = X(j,isfinite(X(j,:)));
        if numel(v) >= 2, xr(j) = max(v)-min(v); end
    end
    dx_trk = max(xr);

    fprintf('\n  x-range in RAW DETECTIONS : %7.1f px = %5.2f cm\n', ...
        dx_det, dx_det*base.mmPerPx/10);
    fprintf('  x-range in TRACKS (best)  : %7.1f px = %5.2f cm\n', ...
        dx_trk, dx_trk*base.mmPerPx/10);
    fprintf('  motion retained by tracking: %.1f%%\n', 100*dx_trk/dx_det);
    if dx_det > 50 && dx_trk < 0.3*dx_det
        fprintf('  >>> DETECTIONS CONTAIN MOTION THAT TRACKS DO NOT.\n');
        fprintf('      Failure is in ASSOCIATION, not detection.\n');
    elseif dx_det <= 50
        fprintf('  >>> Detections themselves are nearly stationary: DETECTION issue.\n');
    else
        fprintf('  >>> Tracks retain the motion: Stage A looks OK for this trial.\n');
    end

    % ---- where does association drop out? -------------------------------
    fprintf('\n  frames with >=8 detections : %d\n', sum(detN>=8));
    fprintf('  frames with >=8 tracked    : %d\n', sum(trkN>=8));
    trkPad = zeros(nFdet,1);
    nCopy  = min(nFtr, nFdet - fvf + 1);
    if nCopy > 0
        trkPad(fvf:fvf+nCopy-1) = trkN(1:nCopy);
    end
    fprintf('  frames with >=1 detection but 0 tracked: %d\n', ...
        sum(detN > 0 & trkPad == 0));

    % ---- plot -----------------------------------------------------------
    figure('Color','w','Name',char(tg),'Position',[60 60 950 800]);

    subplot(3,1,1);
    plot(1:nFdet, detN, '-'); hold on; grid on;
    plot((1:nFtr)+fvf-1, trkN, '-','LineWidth',1.2);
    yline(8,'k:','8 markers');
    legend({'raw detections/frame','tracked markers/frame'},'Location','best');
    ylabel('count'); xlabel('video frame');
    title(sprintf('%s — detection vs tracking yield', tg),'Interpreter','none');

    subplot(3,1,2);
    plot(1:nFdet, detXmin,'.','MarkerSize',4); hold on; grid on;
    plot(1:nFdet, detXmax,'.','MarkerSize',4);
    legend({'min x of detections','max x of detections'},'Location','best');
    ylabel('x (px)'); xlabel('video frame');
    title(sprintf('RAW detection extent — spans %.0f px = %.2f cm', ...
        dx_det, dx_det*base.mmPerPx/10));

    subplot(3,1,3);
    plot((1:nFtr)+fvf-1, X', '-'); grid on;
    ylabel('tracked x (px)'); xlabel('video frame');
    title(sprintf('TRACKED marker x — best marker spans %.0f px = %.2f cm', ...
        dx_trk, dx_trk*base.mmPerPx/10));
end

fprintf('\nRead-only — nothing written.\n');
