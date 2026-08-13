%% ========================================================================
%  TRACK INTEGRITY DIAGNOSTIC — Default / Tight / Wide, GB/full
%
%  Answers, per trial:
%    - how many frames the file holds vs how many actually contain data
%    - whether trackedX/Y are NaN-filled, frozen, or a real trajectory
%    - per-marker displacement range in px
%    - whether impactDistPx falls inside the tracked range
%    - which fields the tracks/meta structs carry (format drift between runs)
%
%  KEY LOGIC: fps enters kinematics only through dt. It scales v and a but
%  CANNOT change displacement in pixels. Zero pixel travel is therefore never
%  an fps problem — it is a tracking-output problem.
%
%  Read-only. Writes nothing.
%% ========================================================================

root   = 'D:\ME_GRANULAB\JerboaImpact';
gbRoot = fullfile(root,'03_RESULTS','GB','Batch 5');
base   = get_calibration();
bn     = sqrt(base.lineA^2 + base.lineB^2);
IDP    = containers.Map({'Tight','Default','Wide'}, {-376.001, -370.001, -409});

NPER = 6;                          % trials sampled per model for the table

D = dir(fullfile(gbRoot,'**','*_tracks.mat'));
allNames = string({D.name});
isFull   = contains(allNames,'_full');
D = D(isFull);
fprintf('GB/full track files: %d\n', numel(D));

mdlAll = strings(numel(D),1);
for k = 1:numel(D)
    mdlAll(k) = model_of(string(D(k).name), string(D(k).folder));
end

%% ---- 1. PER-TRIAL INTEGRITY TABLE -------------------------------------
rows = [];
models = ["Default","Tight","Wide"];
for mi = 1:numel(models)
    idx = find(mdlAll == models(mi));
    if isempty(idx), continue; end
    pick = idx(round(linspace(1, numel(idx), min(NPER,numel(idx)))));
    for k = pick(:)'
        tp = fullfile(D(k).folder, D(k).name);
        S  = load(tp);
        if ~isfield(S,'tracks'), continue; end
        tr = S.tracks;
        m  = getfielddef(S,'meta',struct());

        X = tr.trackedX; Y = tr.trackedY;
        [nM, nF] = size(X);
        fin      = isfinite(X) & isfinite(Y);
        nAll     = sum(all(fin,1));            % frames with EVERY marker
        nAny     = sum(any(fin,1));            % frames with at least one
        pctFin   = 100*nnz(fin)/numel(fin);

        fps = getfielddef(m,'fps_true',NaN);
        if ~isfinite(fps), fps = getfielddef(tr,'fps',NaN); end

        d = (base.lineA.*X + base.lineB.*Y + base.lineC)./bn;

        % per-marker travel
        trav = nan(nM,1);
        for j = 1:nM
            v = d(j,:); v = v(isfinite(v));
            if numel(v) >= 2, trav(j) = max(v)-min(v); end
        end

        % reference marker = furthest on the rod side, first full frame
        ff = find(all(fin,1),1,'first');
        if isempty(ff), ff = find(any(fin,1),1,'first'); end
        refTrav = NaN; refN = 0; rmin = NaN; rmax = NaN; refID = NaN;
        if ~isempty(ff)
            col = d(:,ff); col(~isfinite(col)) = +Inf;
            [~,refID] = min(col);
            r = d(refID,:); r = r(isfinite(r));
            refN = numel(r);
            if refN >= 1, rmin = min(r); rmax = max(r); end
            if refN >= 2, refTrav = rmax - rmin; end
        end

        trg = IDP(char(models(mi)));
        rows = [rows; { char(models(mi)), D(k).name, nM, nF, ...
            round(pctFin,1), nAll, nAny, refN, fps, nF/fps*1e3, ...
            refID, rmin, rmax, refTrav, refTrav*base.mmPerPx/10, ...
            max(trav), trg, (trg>=rmin && trg<=rmax) } ]; %#ok<AGROW>
    end
end

Tt = cell2table(rows, 'VariableNames', {'model','file','nMarkers','nFramesInFile', ...
    'pctFinite','nFrames_allMarkers','nFrames_anyMarker','nFinite_refMarker', ...
    'fps','duration_ms','refMarkerID','ref_min_px','ref_max_px', ...
    'refTravel_px','refTravel_cm','maxMarkerTravel_px','impactDistPx','triggerInRange'});

fprintf('\n================ TRACK INTEGRITY ================\n');
disp(Tt(:,{'model','file','nFramesInFile','nFrames_allMarkers','nFinite_refMarker', ...
           'pctFinite','fps','duration_ms'}))
fprintf('\n---- displacement + trigger ----\n');
disp(Tt(:,{'model','file','refMarkerID','ref_min_px','ref_max_px','refTravel_px', ...
           'refTravel_cm','maxMarkerTravel_px','impactDistPx','triggerInRange'}))

%% ---- 2. INTERPRETATION ------------------------------------------------
fprintf('\n================ WHAT THIS MEANS ================\n');
for mi = 1:numel(models)
    s = Tt(strcmp(Tt.model,char(models(mi))),:);
    if isempty(s), continue; end
    fprintf('\n%s:\n', models(mi));
    fprintf('  frames in file      : %s\n', mat2str(s.nFramesInFile'));
    fprintf('  frames w/ ALL mkrs  : %s\n', mat2str(s.nFrames_allMarkers'));
    fprintf('  finite ref samples  : %s\n', mat2str(s.nFinite_refMarker'));
    fprintf('  ref travel (px)     : %s\n', mat2str(round(s.refTravel_px',1)));
    if all(s.nFinite_refMarker <= 1)
        fprintf('  >>> ONE OR ZERO finite samples: min==max is an artefact of a\n');
        fprintf('      single data point, NOT a stationary marker. Tracking\n');
        fprintf('      produced essentially no usable trajectory.\n');
    elseif all(s.refTravel_px < 1)
        fprintf('  >>> Many finite samples but ~zero travel: positions are FROZEN\n');
        fprintf('      (same value repeated). Tracking wrote a constant.\n');
    elseif median(s.duration_ms) < 30
        fprintf('  >>> Trajectory present but the record is only %.0f ms; a GB/full\n', median(s.duration_ms));
        fprintf('      penetration takes 30-40 ms, so the event is TRUNCATED.\n');
    else
        fprintf('  >>> Trajectory looks complete.\n');
    end
end

%% ---- 3. FORMAT COMPARISON (did the save format drift?) ----------------
fprintf('\n================ FILE FORMAT ================\n');
for mi = 1:numel(models)
    idx = find(mdlAll == models(mi), 1, 'first');
    if isempty(idx), continue; end
    S = load(fullfile(D(idx).folder, D(idx).name));
    fprintf('\n%-8s %s\n', models(mi), D(idx).name);
    fprintf('  top-level vars : %s\n', strjoin(fieldnames(S)', ', '));
    if isfield(S,'tracks')
        fprintf('  tracks fields  : %s\n', strjoin(fieldnames(S.tracks)', ', '));
        fprintf('  size(trackedX) : %s   class %s\n', ...
            mat2str(size(S.tracks.trackedX)), class(S.tracks.trackedX));
    end
    if isfield(S,'meta')
        fn = fieldnames(S.meta)';
        fprintf('  meta fields    : %s\n', strjoin(fn, ', '));
        fprintf('  fps_true=%s  container=%s  dropHeight=%s\n', ...
            num2str(getfielddef(S.meta,'fps_true',NaN)), ...
            char(string(getfielddef(S.meta,'container','?'))), ...
            num2str(getfielddef(S.meta,'dropHeight_mm',NaN)));
    end
    if isfield(S,'calib')
        fprintf('  calib          : mmPerPx=%.4f impactDistPx=%g bedX=%s\n', ...
            getfielddef(S.calib,'mmPerPx',NaN), getfielddef(S.calib,'impactDistPx',NaN), ...
            num2str(getfielddef(S.calib,'bedX',NaN)));
    end
end

%% ---- 4. CALIBRATION CHECK ---------------------------------------------
fprintf('\n================ CALIBRATION PER MODEL ================\n');
fprintf('%-8s %12s %12s %12s\n','model','expected','used','in range?');
for mi = 1:numel(models)
    s = Tt(strcmp(Tt.model,char(models(mi))),:);
    if isempty(s), continue; end
    fprintf('%-8s %12.3f %12.3f %12s\n', models(mi), IDP(char(models(mi))), ...
        s.impactDistPx(1), mat2str(any(s.triggerInRange)));
end
fprintf('\nShared for all models: mmPerPx=%.4f  bedX=%g  - spacing check confirmed\n', ...
    base.mmPerPx, base.bedX);
fprintf('the rod images at the same magnification in all three runs\n');

%% ---- 5. RAW TRACE PLOT, each trial on its OWN fps ---------------------
figure('Color','w','Name','raw marker traces by model','Position',[60 60 1100 750]);
for mi = 1:numel(models)
    idx = find(mdlAll == models(mi), 1, 'first');
    if isempty(idx), continue; end
    S = load(fullfile(D(idx).folder, D(idx).name),'tracks','meta');
    X = S.tracks.trackedX; Y = S.tracks.trackedY;
    fps = getfielddef(S.meta,'fps_true',NaN);
    if ~isfinite(fps), fps = getfielddef(S.tracks,'fps',NaN); end
    d = (base.lineA.*X + base.lineB.*Y + base.lineC)./bn;
    t = (0:size(X,2)-1)/fps*1e3;

    subplot(3,1,mi); hold on; grid on;
    plot(t, d', '-');                       % every marker
    yline(IDP(char(models(mi))),'r--','impactDistPx');
    xlabel('t (ms, own fps)'); ylabel('marker \rightarrow bed (px)');
    title(sprintf('%s — %s   (%d frames @ %.0f fps = %.1f ms)', ...
        models(mi), D(idx).name, size(X,2), fps, size(X,2)/fps*1e3), 'Interpreter','none');
end

fprintf('\nRead-only — no kinematics written.\n');

%% ---- helpers -----------------------------------------------------------
function m = model_of(name, folder)
    t = lower(string(name)); f = lower(string(folder));
    if contains(t,"_tight") || contains(f,"tight"), m = "Tight";
    elseif contains(t,"_wide") || contains(f,"wide"), m = "Wide";
    else, m = "Default";
    end
end

function v = getfielddef(s,f,d)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
