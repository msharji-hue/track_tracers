function P = track_tracers_2(mode, target, opts)
% TRACK_TRACERS_2  Kinematics pass over saved tracking outputs (Stage B).
%
%   Reads the tracking-only results produced by process_trial (the *_tracks.mat
%   files under 03_RESULTS), runs the velocity-first kd_kinematics on each, and
%   writes per-trial kinematics + a crash-safe batch log. It does NOT re-track,
%   read frames, or touch trackedX/trackedY — tracking is already done and
%   FPS-repaired.
%
%   Follows process_trial's naming/output conventions exactly and shares its
%   helpers (resolve_output_root, write_dryrun_report, batch_log_*).
%
%   USAGE
%     track_tracers_2(mode, target, opts)
%       mode   : 'batch' (default) | 'single' | 'rerun'
%       target : batch  -> output root (parent of 03_RESULTS), or '' to prompt
%                single -> a trialTag (e.g. '285mm_T08_dense') or a _tracks.mat path
%                rerun  -> a trialTag; forces overwrite for that trial
%       opts   : .root    output root (if not given as target)
%                .policy  'reuse'(default) | 'resume' | 'overwrite'
%                         use 'overwrite' to FORCE reprocessing of trials
%                         that already have a _kin.mat (otherwise SKIP(done))
%                .limit   process only first N trials (test runs)     [0 = all]
%                .select  cellstr of trialTags to process (test subset)
%                .dryRun  true -> list trials + planned paths, write nothing.
%                         Does NOT compute kinematics; use .preview for that.
%                .preview true -> run the FULL kinematics on every selected
%                         trial and show the figures, but write nothing to
%                         disk: no _kin.mat, no _kin_scalars.csv, no batch log.
%                         Forces .figures to 'show' unless it was set to
%                         'none'. Use this to inspect the whole batch before
%                         committing. Returns the results in the workspace via
%                         the optional output argument.
%                .figures 'none'(batch default) | 'show' | 'save'
%                .massG   carriage+foot mass in grams (default 65) for f in N
%                .model   model-level subfolder under 03_RESULTS to search,
%                         e.g. 'Default model'  ['' = search all models]
%                .saveEventFrames true (default) -> after each trial is
%                         committed, re-export the frames around the impact
%                         ([impact-pad, stop+pad]) from the raw video into the
%                         trial's 01_FRAMES mirror, raw-indexed, so the impact
%                         QA and diag_impact_frame have PNGs to read. Stage A
%                         streams and keeps none. Skipped under dryRun/preview,
%                         idempotent, and never fails a trial.
%                .eventFramePad  [pre post] frames, default [20 20]
%                .rawRoot  capture tree for the re-export
%                         (default $JERBOA_RAW_ROOT, else the campaign-1 tree)
%                .impactCheck true (default) -> after each trial is committed,
%                         write an impact-QA PNG via diag_impact_frame into
%                         03_RESULTS/_batch_logs/impact_checks. Skipped under
%                         dryRun and preview. A QA failure logs a WARN line and
%                         never fails the trial.
%
%   OUTPUTS (per trial, sibling of the tracks/ folder — stays in the SAME tree)
%     <...>/<container>/kinematics/<trialTag>_kin.mat          meta + kin + calib
%     <...>/<container>/kinematics/<trialTag>_kin_scalars.csv  one-row summary
%   Batch log -> 03_RESULTS/[<model>/]_batch_logs/kin_log_<stamp>.txt + kin_progress_<stamp>.csv

    P = table();          % defined up front so every early return is safe
    if nargin < 1 || isempty(mode), mode = 'batch'; end
    if nargin < 2, target = ''; end
    if nargin < 3 || isempty(opts), opts = struct(); end
    opts.root    = getfld(opts,'root','');
    opts.policy  = lower(getfld(opts,'policy','reuse'));
    opts.limit   = getfld(opts,'limit',0);
    opts.select  = getfld(opts,'select',{});
    opts.dryRun  = getfld(opts,'dryRun',false);
    opts.preview = getfld(opts,'preview',false);
    opts.figures = lower(getfld(opts,'figures','none'));
    if opts.preview && strcmp(opts.figures,'none') && ~isfield(opts,'figures')
        opts.figures = 'show';
    end
    opts.massG   = getfld(opts,'massG',65);
    opts.model   = getfld(opts,'model','');   % e.g. 'Default model'; '' = search all
    % Per-trial impact QA. Writes one PNG per committed trial via
    % diag_impact_frame so a batch can be checked by eye afterwards. Never runs
    % under dryRun or preview (nothing is committed in either), and a QA failure
    % never fails the trial -- the kinematics are already on disk by then.
    opts.impactCheck = getfld(opts,'impactCheck',true);
    % Event-frame subset. Stage A no longer keeps 01_FRAMES, so the frames
    % around the impact are re-exported here -- a few dozen per trial rather
    % than the whole clip -- into the standard 01_FRAMES mirror, using the raw
    % frame_%05d.png naming so diag_impact_frame's PNG-first lookup finds them
    % with no change. Skipped under dryRun and preview.
    opts.saveEventFrames = getfld(opts,'saveEventFrames',true);
    opts.eventFramePad   = getfld(opts,'eventFramePad',[20 20]);
    opts.eventFramePad   = double(opts.eventFramePad(:).');
    if numel(opts.eventFramePad) == 1
        opts.eventFramePad = opts.eventFramePad([1 1]);
    end
    opts.rawRoot = getfld(opts,'rawRoot', local_default_raw_root());
    mode = lower(strtrim(mode));

    thisDir = fileparts(mfilename('fullpath'));
    codeDir = fileparts(thisDir);
    addpath(fullfile(codeDir,'src'));
    calib = get_calibration();

    % ── resolve output root ───────────────────────────────────────────────
    if isempty(opts.root) && strcmp(mode,'batch') && ~isempty(target) ...
            && ~endsWith(lower(target),'.mat')
        opts.root = target;
    end
    root = resolve_output_root(opts.root);
    if isempty(root), fprintf('Cancelled.\n'); return; end
    resultsRoot = fullfile(root,'03_RESULTS');
    if ~isempty(opts.model)
        resultsRoot = fullfile(resultsRoot, opts.model);   % .../03_RESULTS/Default model
    end
    if ~isfolder(resultsRoot)
        error('Results folder not found:\n  %s\nCheck ROOT and the ''model'' option.', resultsRoot);
    end
    logDir = fullfile(resultsRoot,'_batch_logs');           % logs stay inside the model tree

    % ── discover trials ───────────────────────────────────────────────────
    switch mode
        case 'batch'
            items     = discover_tracks(resultsRoot);
            inputDesc = resultsRoot;
        case {'single','rerun'}
            items     = resolve_single(target, resultsRoot);
            inputDesc = target;
            if strcmp(mode,'rerun'), opts.policy = 'overwrite'; end
        otherwise
            error('Unknown mode "%s" (batch|single|rerun)', mode);
    end

    % ── filter (select / limit) ───────────────────────────────────────────
    if ~isempty(opts.select) && ~isempty(items)
        items = items(ismember({items.trialTag}, opts.select));
    end
    if opts.limit > 0 && numel(items) > opts.limit
        items = items(1:opts.limit);
    end
    if isempty(items), warning('No trials matched.'); return; end

    % ── dry run ───────────────────────────────────────────────────────────
    if opts.dryRun
        rows = repmat(struct('head','','ok',false,'pathLines',{{}}), numel(items), 1);
        for i = 1:numel(items)
            it  = items(i);
            act = 'PROC';
            if it.exists && any(strcmp(opts.policy,{'reuse','resume'})), act = 'SKIP(done)'; end
            if it.exists && strcmp(opts.policy,'overwrite'),             act = 'OVERWRITE';  end
            rows(i).head = sprintf('[%3d] %-10s %-22s  cond=%s drop=%gmm fps=%s', ...
                i, act, it.trialTag, it.condition, it.dropHeight_mm, num2str(it.fps,'%.1f'));
            rows(i).ok        = true;
            rows(i).pathLines = {'tracks', it.tracksPath; 'kin', it.kinMat};
        end
        write_dryrun_report(logDir, inputDesc, rows);
        P = table();
        return;
    end

    PREVIEW = struct('trialTag',{},'condition',{},'dropHeight_mm',{}, ...
                     'fps',{},'impactDistPx',{},'impact_index',{}, ...
                     'stopFrame',{},'v0_cm_s',{},'d_final_cm',{},'t_stop_s',{});
    P = table();

    % ── process ───────────────────────────────────────────────────────────
    if opts.preview
        fprintf(['\nPREVIEW: computing kinematics for %d trial(s). Nothing will be\n' ...
                 'written -- no _kin.mat, no scalars CSV, no batch log.\n\n'], numel(items));
    end
    stamp  = datestr(now,'yyyymmdd_HHMMSS');
    header = sprintf([ ...
        '============================================================\n' ...
        ' KINEMATICS LOG  %s\n Input  : %s\n' ...
        ' Policy : %s   Limit : %d   Figures : %s   mass : %g g\n' ...
        ' Total trials : %d\n' ...
        '============================================================\n\n'], ...
        stamp, inputDesc, opts.policy, opts.limit, opts.figures, opts.massG, numel(items));
    csvHeader = sprintf(['idx,trialTag,condition,dropHeight_mm,fps,status,', ...
        'v0_cm_s,d_final_cm,t_stop_s,a_stop_cm_s2,reason,kinPath\n']);
    L = batch_log_init(logDir, stamp, {'kin_log_','kin_progress_','kin_retry_'}, header, csvHeader);

    nOK = 0; nSkip = 0; nFail = 0; nQA = 0; nQAfail = 0; nEF = 0; nEFfail = 0; total = numel(items);
    for i = 1:total
        it = items(i);
        fprintf('[%3d/%3d] %s ...\n', i, total, it.trialTag);

        if it.exists && any(strcmp(opts.policy,{'reuse','resume'}))
            nSkip = nSkip + 1;
            log_kin_row(L, i, total, it, 'SKIPPED', 'already has _kin.mat', []);
            continue;
        end

        try
            % 'calib' is requested too: Stage A saves the per-model
            % calibration alongside the tracks, and it is used below.
            S = load(it.tracksPath, 'meta', 'tracks', 'calib');
        catch ME
            nFail = nFail + 1;
            log_kin_row(L, i, total, it, 'FAILED', ['load: ' ME.message], []);
            continue;
        end
        meta = S.meta; tracks = S.tracks;

        fps = resolve_fps(it.tracksPath, meta, tracks);
        if ~isfinite(fps)
            nFail = nFail + 1;
            log_kin_row(L, i, total, it, 'FAILED', ...
                'no plausible fps in [200,20000] from scalars CSV / meta / tracks', []);
            continue;
        end
        it.fps = fps;
        it.calibSource = 'global';
        if isfield(S,'calib') && isstruct(S.calib) && isfield(S.calib,'impactDistPx')
            it.calibSource   = 'tracks';
            it.impactDistPx  = S.calib.impactDistPx;
        else
            it.impactDistPx  = NaN;
        end

        try
            % Calibration comes from the TRACKS FILE when it is there.
            % Stage A saves the exact calib it used, so whatever bed line and
            % search-window centre produced the tracks flows through
            % automatically and there is no second place to keep in sync.
            % Trials predating the model work have no saved calib and fall back
            % to the global one -- which is the value they were tracked with
            % anyway.
            %
            % From the 2026-08 campaign impactDistPx is a standardized -360 px
            % for every model, so the fallback no longer risks applying one
            % model's trigger to another. Older tracks files still carry their
            % original per-model value, and reading it from the file keeps those
            % trials reproducible.
            if isfield(S,'calib') && isstruct(S.calib) && ...
               isfield(S.calib,'impactDistPx')
                calibT = S.calib;
            else
                calibT = get_calibration();
            end
            kin = kd_kinematics(tracks.trackedX, tracks.trackedY, calibT, 1/fps);

            % Rod rotation diagnostic: line through marker centres vs the
            % reference (first full) frame. Scalars over impact:stop; noise
            % floor from the quiet post-stop frames.
            try
                [angDeg, angS] = rod_angle(tracks.trackedX, tracks.trackedY, ...
                    'impact', kin.impact_index, 'stop', kin.stopFrame);
                angS.delta_deg = angDeg(:);
                kin.rodAngle   = angS;
            catch
                kin.rodAngle = struct('peak_abs_deg',NaN,'peak_signed_deg',NaN, ...
                    'net_deg',NaN,'range_deg',NaN,'noise_sd_deg',NaN, ...
                    'peakFrame',NaN,'refFrame',NaN,'delta_deg',[]);
            end
        catch ME
            nFail = nFail + 1;
            log_kin_row(L, i, total, it, 'FAILED', ['kd_kinematics: ' ME.message], []);
            continue;
        end

        fprintf('    calib: %s  impactDistPx = %g  (fps %.0f)\n', ...
                it.calibSource, getfld(calibT,'impactDistPx',NaN), fps);


        calib = calibT;                          % the calibration actually used
        if opts.preview
            % Everything above has run; nothing below touches the disk.
            PREVIEW(end+1) = struct('trialTag',it.trialTag, ...
                'condition',it.condition, 'dropHeight_mm',it.dropHeight_mm, ...
                'fps',fps, 'impactDistPx',getfld(calib,'impactDistPx',NaN), ...
                'impact_index',getfld(kin,'impact_index',NaN), ...
                'stopFrame',getfld(kin,'stopFrame',NaN), ...
                'v0_cm_s',getfld(kin,'v0_cm_s',NaN), ...
                'd_final_cm',getfld(kin,'d_final_cm',NaN), ...
                't_stop_s',getfld(kin,'t_stop_s',NaN)); %#ok<AGROW>
        else
            if ~exist(it.kinDir,'dir'), mkdir(it.kinDir); end

            % ── event-frame subset ───────────────────────────────────────
            % Before the scalars are written, so eventFrameRange can be
            % recorded in both the .mat and the CSV. Same windowStart +
            % firstValidFrame arithmetic diag_impact_frame uses.
            meta.eventFrameRange = [NaN NaN];
            if opts.saveEventFrames
                try
                    meta.eventFrameRange = export_event_frames(it, meta, kin, opts);
                    nEF = nEF + 1;
                catch ME_EF
                    nEFfail = nEFfail + 1;
                    warnLine = sprintf(['  WARN event frames skipped for %s: %s\n' ...
                                        '       (kinematics were written; frames only)\n'], ...
                                       it.trialTag, ME_EF.message);
                    fprintf('%s', warnLine);
                    batch_log_row(L, warnLine, '');
                end
            end

            save(it.kinMat, 'meta', 'kin', 'calib', '-v7.3');
            write_kin_scalars(it, meta, kin, calib, opts.massG);

            % ── per-trial impact QA ──────────────────────────────────────
            % Strictly after the _kin.mat exists: diag_impact_frame reads it
            % back off disk. Wrapped so a QA problem -- a purged frame, a
            % missing raw clip, a display fault -- can never cost a trial whose
            % kinematics are already committed.
            if opts.impactCheck
                try
                    diag_impact_frame(it.trialTag, 'Root', root, ...
                                      'Save', true, 'Show', false);
                    nQA = nQA + 1;
                catch ME_QA
                    nQAfail = nQAfail + 1;
                    warnLine = sprintf(['  WARN impact QA skipped for %s: %s\n' ...
                                        '       (kinematics were written; QA only)\n'], ...
                                       it.trialTag, ME_QA.message);
                    fprintf('%s', warnLine);
                    batch_log_row(L, warnLine, '');
                end
            end
        end

        if ~strcmp(opts.figures,'none')
            f1 = fig_triptych(kin, meta);
            if strcmp(opts.figures,'save')
                saveas(f1, fullfile(it.kinDir,[it.trialTag '_kin_triptych.png']));
                close(f1);
            end
        end

        nOK = nOK + 1;
        log_kin_row(L, i, total, it, 'OK', '', kin);
    end

    fin = sprintf([ ...
        '============================================================\n' ...
        ' DONE  total=%d  OK=%d  SKIPPED=%d  FAILED=%d\n' ...
        '============================================================\n'], ...
        total, nOK, nSkip, nFail);
    batch_log_finalize(L, fin);
    fprintf('\nKinematics batch done: OK=%d  SKIPPED=%d  FAILED=%d\nLog: %s\n', ...
        nOK, nSkip, nFail, L.txt);

    if opts.saveEventFrames && ~opts.preview
        fprintf('Event frames: %d trial(s) exported, %d skipped\n', nEF, nEFfail);
    end
    if opts.impactCheck && ~opts.preview
        % diag_impact_frame writes under <root>/03_RESULTS/_batch_logs, with no
        % model level, so this is deliberately built from root and not from
        % resultsRoot (which carries the model subfolder when opts.model is set).
        fprintf('Impact QA: %d written, %d skipped\n  %s\n', nQA, nQAfail, ...
                fullfile(root,'03_RESULTS','_batch_logs','impact_checks'));
    end

    if opts.preview && ~isempty(PREVIEW)
        P = struct2table(PREVIEW);
        fprintf('\n=== PREVIEW SUMMARY (nothing written) ===\n');
        fprintf('  %d trial(s) computed\n', height(P));
        fprintf('  v0      %.1f - %.1f cm/s\n', min(P.v0_cm_s), max(P.v0_cm_s));
        fprintf('  d_final %.3f - %.3f cm\n', min(P.d_final_cm), max(P.d_final_cm));
        fprintf('  fps     %.0f - %.0f\n', min(P.fps), max(P.fps));
        u = unique(P.impactDistPx);
        fprintf('  impactDistPx in use: %s\n', strjoin(compose('%g',u'), ', '));
        fprintf(['\n  Inspect P, then re-run without ''preview'' to write:\n' ...
                 '    track_tracers_2(''batch'', root, struct(''figures'',''none''))\n\n']);
    end
end

% ═════════════════════════════ helpers ═══════════════════════════════════
function v = getfld(s, f, dflt)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
function v = fld(s, f)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end

function it = build_item(tp, m)
    [~, base] = fileparts(tp);
    tag = getfld(m,'trialTag', strrep(base,'_tracks',''));
    it.trialTag      = tag;
    it.material      = getfld(m,'material','');
    it.container     = getfld(m,'container','');
    it.condition     = sprintf('%s/%s', it.material, it.container);
    it.dropHeight_mm = getfld(m,'dropHeight_mm',NaN);
    it.trialNum      = getfld(m,'trialNum',NaN);
    it.fps           = resolve_fps(tp, m, []);
    it.tracksPath    = tp;
    it.resultsDir    = fileparts(fileparts(tp));           % .../<container>
    it.kinDir        = fullfile(it.resultsDir,'kinematics');
    it.kinMat        = fullfile(it.kinDir,[tag '_kin.mat']);
    it.exists        = exist(it.kinMat,'file') > 0;
end

function items = discover_tracks(resultsRoot)
    D = dir(fullfile(resultsRoot,'**','*_tracks.mat'));
    items = repmat(build_item('seed', struct()), 0, 1);   % 0x1 with right fields
    for k = 1:numel(D)
        tp = fullfile(D(k).folder, D(k).name);
        try, Sm = load(tp,'meta'); m = Sm.meta; catch, continue; end
        items(end+1) = build_item(tp, m); %#ok<AGROW>
    end
    if ~isempty(items)
        keys = arrayfun(@(x) sprintf('%s_%08.1f_%04.0f', x.condition, ...
            x.dropHeight_mm, x.trialNum), items, 'UniformOutput', false);
        [~, ord] = sort(keys);
        items = items(ord);
    end
end

function items = resolve_single(target, resultsRoot)
    if endsWith(lower(target),'.mat') && exist(target,'file')
        tp = target;
    else
        D = dir(fullfile(resultsRoot,'**',[target '_tracks.mat']));
        if isempty(D)
            error('No _tracks.mat found for "%s" under %s', target, resultsRoot);
        end
        tp = fullfile(D(1).folder, D(1).name);
    end
    Sm = load(tp,'meta');
    items = build_item(tp, Sm.meta);
end

function write_kin_scalars(it, meta, kin, calib, massG) %#ok<INUSD>
    p  = fullfile(it.kinDir,[it.trialTag '_kin_scalars.csv']);
    fid = fopen(p,'w');
    ra = getfld(kin,'rodAngle',struct());
    fprintf(fid, ['material,batch,dropHeight_mm,', ...
        'trialNum,condition,phi,fps,', ...
        'impact_frame,stop_frame,v0_cm_s,d_final_cm,a_stop_cm_s2,t_stop_s,', ...
        'impactDistPx,bedX,windowStart,windowEnd,autoWindow,', ...
        'eventFrameFirst,eventFrameLast,', ...
        'rodAngle_peak_abs_deg,rodAngle_peak_signed_deg,rodAngle_net_deg,', ...
        'rodAngle_range_deg,rodAngle_noise_sd_deg,rodAngle_refFrame\n']);
    fprintf(fid, ['%s,%s,%g,%g,%s,%.4f,%.4f,%d,%d,%.4f,%.4f,%.2f,%.6f,%g,%g,', ...
        '%g,%g,%d,%g,%g,', ...
        '%.4f,%.4f,%.4f,%.4f,%.4f,%g\n'], ...
        getfld(meta,'material',''), getfld(meta,'batchName',''), ...
        getfld(meta,'dropHeight_mm',NaN), getfld(meta,'trialNum',NaN), it.condition, ...
        fld(meta,'phi'), it.fps, kin.impact_index, kin.stopFrame, ...
        kin.v0_cm_s, kin.d_final_cm, getfld(kin,'a_stop_cm_s2',NaN), kin.t_stop_s, ...
        getfld(calib,'impactDistPx',NaN), getfld(calib,'bedX',NaN), ...
        getfld(meta,'windowStart',NaN), getfld(meta,'windowEnd',NaN), ...
        double(isequal(getfld(meta,'autoWindow',false), true)), ...
        local_ef(meta,1), local_ef(meta,2), ...
        getfld(ra,'peak_abs_deg',NaN), getfld(ra,'peak_signed_deg',NaN), ...
        getfld(ra,'net_deg',NaN), getfld(ra,'range_deg',NaN), ...
        getfld(ra,'noise_sd_deg',NaN), getfld(ra,'refFrame',NaN));
    fclose(fid);
end

function log_kin_row(L, i, total, it, status, reason, kin)
    if isempty(kin)
        v0=NaN; dd=NaN; ts=NaN; as_=NaN; kp='';
    else
        v0=kin.v0_cm_s; dd=kin.d_final_cm; ts=kin.t_stop_s;
        as_=getfld(kin,'a_stop_cm_s2',NaN); kp=it.kinMat;
    end
    txt = sprintf('[%3d/%3d] %-22s  %s\n', i, total, it.trialTag, status);
    txt = [txt sprintf('          cond=%s  drop=%gmm  fps=%s\n', ...
        it.condition, it.dropHeight_mm, num2str(it.fps,'%.1f'))];
    if ~isempty(reason), txt = [txt sprintf('          note   : %s\n', reason)]; end
    if strcmp(status,'OK')
        txt = [txt sprintf('          v0=%.1f cm/s  d=%.3f cm  a_stop=%.1f\n', v0,dd,as_)];
        txt = [txt sprintf('          kin    : %s\n', kp)];
    end
    txt = [txt sprintf('\n')];

    csvRow = sprintf('%d,%s,%s,%g,%s,%s,%.4f,%.4f,%.6f,%.2f,"%s",%s\n', ...
        i, it.trialTag, it.condition, it.dropHeight_mm, num2str(it.fps,'%.1f'), ...
        status, v0, dd, ts, as_, reason, kp);

    retryLine = '';
    if strcmp(status,'FAILED'), retryLine = it.tracksPath; end
    batch_log_row(L, txt, csvRow, retryLine);
end

% ── figures (test/single only) ────────────────────────────────────────────
function f = fig_triptych(kin, meta)
    tag = getfld(meta,'trialTag','?');
    f = figure('Name',['Kinematics — ' tag],'Color','w');
    t = kin.t_s;
    subplot(3,1,1); plot(t, kin.z, '-'); grid on; ylabel('z (cm)');
    title(sprintf('%s   (impact=%d, stop=%d)', tag, kin.impact_index, kin.stopFrame), ...
        'Interpreter','none');
    xline(0,'--'); xline(kin.t_stop_s,'--r');
    subplot(3,1,2); plot(t, kin.v, '-'); grid on; ylabel('v (cm/s)');
    xline(0,'--'); xline(kin.t_stop_s,'--r');
    subplot(3,1,3); plot(t, kin.a_plus_g, '-'); grid on;
    ylabel('a+g (cm/s^2)'); xlabel('t (s)'); yline(0,':');
    xline(0,'--'); xline(kin.t_stop_s,'--r');
end
function rng = export_event_frames(it, meta, kin, opts)
%EXPORT_EVENT_FRAMES  Re-export the frames around the impact for ONE trial.
%
%   Stage A streams and keeps no PNGs, so the handful of frames actually needed
%   for visual QA are written here instead of the whole clip. They go to the
%   trial's standard 01_FRAMES mirror with the SAME raw-indexed frame_%05d.png
%   naming export_frames uses, so diag_impact_frame's PNG-first lookup finds
%   them with no change at all.
%
%   Range, using exactly the arithmetic diag_impact_frame uses:
%       tracking index k  ->  raw frame  windowStart + firstValidFrame + k - 2
%   so
%       [rawImpact - padPre, rawStop + padPost]   clamped to the video.
%
%   Idempotent: if every PNG in the range is already present the export is
%   skipped and the range still returned, so re-running a batch is cheap.

    wStart = 1;
    if isfield(meta,'windowStart') && isfinite(meta.windowStart)
        wStart = meta.windowStart;
    end
    base      = wStart + meta.firstValidFrame - 1;   % raw frame of tracking index 1
    rawImpact = base + kin.impact_index - 1;
    rawStop   = base + kin.stopFrame    - 1;

    lo = rawImpact - opts.eventFramePad(1);
    hi = rawStop   + opts.eventFramePad(2);

    % Frames mirror the results path under 01_FRAMES, honouring the scratch
    % root exactly as process_trial's build_leaf_dirs does.
    resultDir = fileparts(fileparts(it.kinMat));     % <resultsDir>/kinematics -> <resultsDir>
    if ~contains(resultDir, '03_RESULTS')
        error('export_event_frames:noResultsRoot', ...
              'Cannot derive the 01_FRAMES mirror from %s', resultDir);
    end
    framesRoot = getenv('JERBOA_FRAMES_ROOT');
    if isempty(framesRoot)
        framesDir = strrep(resultDir, '03_RESULTS', '01_FRAMES');
    else
        rel       = extractAfter(string(resultDir), "03_RESULTS");
        framesDir = char(fullfile(framesRoot, '01_FRAMES', strip(rel, filesep)));
    end

    % Model is required: the same stem exists under every model folder.
    mdl = '';
    if isfield(meta,'model'), mdl = meta.model; end
    videoPath = find_raw_video(opts.rawRoot, meta, mdl);
    v         = open_video(videoPath);
    nAvail    = floor(v.Duration * v.FrameRate);
    lo        = max(1, round(lo));
    hi        = min(nAvail, round(hi));
    if hi < lo
        error('export_event_frames:emptyRange', ...
              'Event range [%d, %d] is empty (video has %d frames).', lo, hi, nAvail);
    end

    % Already there? Then nothing to do.
    if isfolder(framesDir)
        present = true;
        for f = lo:hi
            if ~isfile(fullfile(framesDir, sprintf('frame_%05d.png', f)))
                present = false; break
            end
        end
        if present
            fprintf('  event frames already present (%d-%d)\n', lo, hi);
            rng = [lo hi];
            return
        end
    end

    % Same filter Stage A used, so these PNGs match what it detected on.
    filt = 'sharpen';
    if isfield(meta,'provenance') && isfield(meta.provenance,'filterType') ...
            && ~isempty(meta.provenance.filterType)
        filt = meta.provenance.filterType;
    end

    if ~isfolder(framesDir), mkdir(framesDir); end
    export_frames(v, lo, hi, framesDir, filt);
    fprintf('  event frames %d-%d -> %s\n', lo, hi, framesDir);
    rng = [lo hi];
end

function v = local_ef(meta, k)
%LOCAL_EF  Element k of meta.eventFrameRange, or NaN when frames were not saved.
    v = NaN;
    if isfield(meta,'eventFrameRange') && numel(meta.eventFrameRange) >= k
        v = meta.eventFrameRange(k);
    end
end

function r = local_default_raw_root()
%LOCAL_DEFAULT_RAW_ROOT  Where the raw capture tree lives.
%   Set JERBOA_RAW_ROOT to point the toolchain at the current campaign without
%   editing code; the fallback is the campaign-1 tree.
    r = getenv('JERBOA_RAW_ROOT');
    if isempty(r), r = 'D:\ME_GRANULAB\Test Batches'; end
end
