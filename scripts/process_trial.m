function process_trial(mode, inputTarget, outputRoot, opts)
% PROCESS_TRIAL  Unified trial-processing entry point.
%
%   Modes:
%     1 - Single trial / manual   GUI pickers + dialogs, debug vars, QA fig
%     2 - Batch                   recursive scan of GB/CHIN x full/shallow
%     3 - Rerun selected          re-track / re-detect / full on one trial
%
%   Interactive (IDE):
%     process_trial                          % menu
%     process_trial('batch')                 % batch with option dialogs
%
%   Headless (terminal):
%     matlab -batch "process_trial('batch','/path/Videos','/path/Out')"
%     matlab -batch "process_trial('batch','/path/Videos','/path/Out', ...
%                     struct('dryRun',true))"
%     matlab -batch "process_trial('batch','/path/Videos','/path/Out', ...
%                     struct('limit',3,'policy','resume','batchLabel','Batch5'))"
%
%   opts (struct, all optional):
%     .dryRun        list videos + metadata + planned paths, process nothing
%     .limit         process only first N videos (0 = all)  [test batch]
%     .batchLabel    inserted into output path; '' = omit the batch level
%     .policy        'reuse'(default) | 'resume' | 'retry' | 'overwrite'
%                      reuse     : reuse frames+detections if present, redo tracks
%                      resume    : additionally SKIP trials already having tracks
%                      retry     : skip trials that already have tracks, and for
%                                  those that do NOT, re-detect from the frames
%                                  instead of reloading the cached detections.
%                                  Use this to retry FAILED trials: 'resume'
%                                  reloads the very detections that failed and
%                                  so reproduces the failure exactly.
%                      overwrite : re-export, re-detect, re-track everything
%     .videoListFile a .txt of full video paths to process (e.g. a retry list)
%     .detectParams  struct of detection-parameter overrides for the FIRST pass
%     .backupParams  struct of detection parameters for a SECOND pass, tried
%                    automatically when the first pass finds no frame with the
%                    expected marker count. Set to struct() to disable. The
%                    default leans the detector toward finding more circles:
%                    higher sensitivity, lower edge threshold, wider radii.
%     .model         foot model: 'Default Model' | 'Tight Model' | 'Wide Model'
%                    When set, the model is appended to trialTag, inserted as a
%                    level in the output path, and used to select the per-model
%                    calibration via get_calibration_model. When '' (default)
%                    every behaviour below is exactly as it was before models
%                    existed, so old runs reproduce bit for bit.

    if nargin < 1, mode        = ''; end
    if nargin < 2, inputTarget = ''; end
    if nargin < 3, outputRoot  = ''; end
    if nargin < 4, opts        = struct(); end
    fromMenu = isempty(mode);

    CFG = struct();
    thisDir        = fileparts(mfilename('fullpath'));
    CFG.codeDir    = fileparts(thisDir);
    CFG.outputRoot = outputRoot;     % '' => prompt
    % Optional: hardcode to skip the prompt:
    % CFG.outputRoot = '/Users/muhannadalsharji/Dropbox/JerboaImpact';
    CFG.nExpectedMarkers = 8;
    CFG.filterType       = 'sharpen';

    addpath(fullfile(CFG.codeDir, 'src'));
    CFG.calib  = get_calibration();
    CFG.params = default_detect_params();

    opts          = normalize_opts(opts);
    opts.fromMenu = fromMenu;

    % Optional detection-parameter override (e.g. Pass 2 recovery). Only the
    % named fields are replaced; everything else in default_detect_params()
    % is preserved, so this cannot silently alter unrelated settings.
    if isfield(opts,'detectParams') && isstruct(opts.detectParams)
        f = fieldnames(opts.detectParams);
        for i = 1:numel(f)
            CFG.params.(f{i}) = opts.detectParams.(f{i});
            fprintf('  [override] params.%s = %g\n', f{i}, opts.detectParams.(f{i}));
        end
    end

    if isempty(mode)
        sel = menu('Select workflow', ...
            '1 - Single trial / manual', ...
            '2 - Batch', ...
            '3 - Rerun selected trial/folder');
        if sel == 0, fprintf('Cancelled.\n'); return; end
        modes = {'single','batch','rerun'};
        mode  = modes{sel};
    end
    mode = normalize_mode(mode);

    switch mode
        case 'single', run_single(CFG, opts);
        case 'batch',  run_batch(CFG, inputTarget, opts);
        case 'rerun',  run_rerun(CFG, inputTarget, opts);
        otherwise,     error('Unknown mode "%s" (single|batch|rerun)', mode);
    end
end

% ═════════════════════════════════════════════════════════════════════════
%  MODE 1 — SINGLE / MANUAL
% ═════════════════════════════════════════════════════════════════════════
function run_single(CFG, opts)
    close all; clc;
    fprintf('\n=== SINGLE TRIAL / MANUAL ===\n');

    [file, inDir] = uigetfile({'*.avi;*.AVI;*.mp4;*.MP4;*.mov;*.MOV','Video Files'}, 'Select a video');
    if isequal(file, 0), fprintf('Cancelled.\n'); return; end
    it = infer_item(fullfile(inDir, file));

    defH = ''; defT = '';
    if ~isnan(it.dropHeight_mm), defH = num2str(it.dropHeight_mm); end
    if ~isnan(it.trialNum),      defT = num2str(it.trialNum);      end
    fields = inputdlg( ...
        {'Material (GB/CHIN)','Batch (blank=none)','Drop height (mm)', ...
         'Trial number','Container (full/shallow)', ...
         'Model (blank=none; Default/Tight/Wide Model)'}, ...
        'Confirm trial info', 1, ...
        {it.material, '', defH, defT, it.container, opts.model});
    if isempty(fields), fprintf('Cancelled.\n'); return; end

    m = build_meta(fields{1}, normalize_batch_label(fields{2}), ...
                   str2double(fields{3}), str2double(fields{4}), ...
                   fields{5}, fields{6});
    outRoot = resolve_output_root(CFG.outputRoot);
    if isempty(outRoot), fprintf('Cancelled.\n'); return; end
    [framesDir, detDir, resultsDir] = build_leaf_dirs(outRoot, m);

    cfg = base_cfg(CFG, m, fullfile(inDir, file), framesDir, detDir, resultsDir);
    cfg.interactive = false;
    cfg.makeQA      = isfield(opts,'makeQA') && isequal(opts.makeQA, true);
    cfg = apply_auto_window(cfg, opts, m.trialTag);
    fprintf('Processing single trial: %s\n', m.trialTag);

    try
        st = process_one_trial(cfg);
        if ~st.ok && st.partial
            fprintf('PARTIAL: %s\nReason: %s (raw tracks saved)\n', file, st.reason);
        elseif ~st.ok
            fprintf('FAILED: %s\nReason: %s\n', file, st.reason);
        end
    catch ME
        fprintf('FAILED: %s\nReason: %s\n', file, ME.message);
    end
end

% ═════════════════════════════════════════════════════════════════════════
%  MODE 2 — BATCH
% ═════════════════════════════════════════════════════════════════════════
function run_batch(CFG, inputRoot, opts)
    clc; fprintf('\n=== BATCH WORKFLOW ===\n');

    % ── Options dialog (interactive only) ────────────────────────────────
    if opts.fromMenu
        a = inputdlg({'Batch label (blank = none in path)','Limit (0 = all)', ...
                      'Model (blank = none)'}, ...
                     'Batch options', 1, {opts.batchLabel, num2str(opts.limit), opts.model});
        if isempty(a), fprintf('Cancelled.\n'); return; end
        opts.batchLabel = normalize_batch_label(a{1});
        opts.limit      = max(0, round(str2double(a{2})));
        opts.model      = strtrim(a{3});

        pol = menu('Redo policy', ...
            'reuse  : keep frames+detections, recompute tracks (default)', ...
            'resume : skip trials already done', ...
            'retry  : skip done trials, RE-DETECT the ones that failed', ...
            'overwrite : redo everything');
        if pol == 0, fprintf('Cancelled.\n'); return; end
        pols = {'reuse','resume','retry','overwrite'}; opts.policy = pols{pol};

        dr = menu('Run type', 'Dry run (list only, no processing)', 'Process for real');
        if dr == 0, fprintf('Cancelled.\n'); return; end
        opts.dryRun = (dr == 1);
    end

    if ~isempty(opts.model)
        fprintf('Model: %s  (tag suffix + path level + per-model calibration)\n', opts.model);
    end

    % ── Gather videos ────────────────────────────────────────────────────
    if ~isempty(opts.videoListFile)
        items     = items_from_list(opts.videoListFile);
        inputDesc = opts.videoListFile;
    else
        if isempty(inputRoot)
            inputRoot = uigetdir(pwd, 'Select ROOT folder (contains GB/ and CHIN/)');
            if isequal(inputRoot, 0), fprintf('Cancelled.\n'); return; end
        end
        items     = scan_video_tree(inputRoot);
        inputDesc = inputRoot;
    end
    total0 = numel(items);
    if total0 == 0, fprintf('No videos found under %s\n', inputDesc); return; end

    if opts.limit > 0 && opts.limit < numel(items)
        items = items(1:opts.limit);
        fprintf('TEST BATCH: limiting to first %d of %d videos.\n', opts.limit, total0);
    end

    outRoot = resolve_output_root(CFG.outputRoot);
    if isempty(outRoot), fprintf('Cancelled.\n'); return; end

    % ── Dry run ──────────────────────────────────────────────────────────
    if opts.dryRun
        print_and_save_dryrun(items, outRoot, opts.batchLabel, inputDesc, opts.model);
        return;
    end

    % ── Process (crash-safe incremental log) ─────────────────────────────
    total = numel(items);
    if opts.autoWindow
        fprintf('autoWindow: ON  (pad %d before / %d after the red-marker span)\n', ...
                opts.windowPad(1), opts.windowPad(2));
    else
        fprintf('autoWindow: OFF (exporting every frame)\n');
    end
    fprintf('Processing %d videos  (policy=%s, batch="%s", model="%s")\n\n', ...
            total, opts.policy, opts.batchLabel, opts.model);
    L = log_init(fullfile(outRoot,'03_RESULTS','_batch_logs'), inputDesc, total, opts);

    nOK = 0; nPartial = 0; nFailed = 0; nSkipped = 0;
    for i = 1:total
        it = items(i);
        fprintf('--------------------------------------------------------\n');
        fprintf('Processing video %d of %d: %s\n', i, total, it.name);

        if ~it.ok
            nFailed = nFailed + 1;
            r = row_from_item(it, opts.batchLabel, 'FAILED', it.reason, NaN, NaN, '');
            fprintf('FAILED: %s\nReason: %s\nContinuing to next video...\n\n', it.name, it.reason);
            log_row(L, i, total, r); continue;
        end

        m = build_meta(it.material, opts.batchLabel, it.dropHeight_mm, it.trialNum, ...
                       it.container, opts.model);
        [framesDir, detDir, resultsDir] = build_leaf_dirs(outRoot, m);

        if any(strcmp(opts.policy,{'resume','retry'})) && ...
           exist(fullfile(resultsDir,'tracks',[m.trialTag '_tracks.mat']),'file')
            nSkipped = nSkipped + 1;
            r = row_from_item(it, opts.batchLabel, 'SKIPPED', 'already processed', NaN, NaN, resultsDir);
            fprintf('SKIPPED (already done): %s\n\n', it.name);
            log_row(L, i, total, r); continue;
        end

        cfg = base_cfg(CFG, m, it.fullpath, framesDir, detDir, resultsDir);
        cfg.interactive = false; cfg.makeQA = true;
        cfg = apply_auto_window(cfg, opts, it.name);
        switch opts.policy
            case 'overwrite'
                cfg.reuseFrames = false; cfg.reuseDetections = false;
            case 'retry'
                % Frames are fine; the DETECTIONS are what failed. Reusing them
                % would reproduce the failure exactly, so force a fresh pass.
                cfg.reuseFrames = true;  cfg.reuseDetections = false;
            otherwise
                cfg.reuseFrames = true;  cfg.reuseDetections = true;
        end

        try
            st = process_one_trial(cfg);

            % ── Pass 2: automatic recovery on a marker-count failure ──────
            % Only fires when pass 1 found no frame with the expected marker
            % count, and only if backup parameters were supplied. Everything
            % else (a read error, a QA failure) is left alone.
            if ~st.ok && ~st.partial && ~isempty(fieldnames(opts.backupParams)) ...
                      && contains(lower(st.reason), 'detected markers')
                fprintf('  Pass 1 found no valid frame. Retrying with backup parameters:\n');
                cfg2 = cfg;
                cfg2.reuseFrames = true; cfg2.reuseDetections = false;
                bf = fieldnames(opts.backupParams);
                for bi = 1:numel(bf)
                    cfg2.params.(bf{bi}) = opts.backupParams.(bf{bi});
                    fprintf('    %s = %s\n', bf{bi}, mat2str(opts.backupParams.(bf{bi})));
                end
                st2 = process_one_trial(cfg2);
                if st2.ok || st2.partial
                    fprintf('  Pass 2 SUCCEEDED.\n');
                    st = st2;
                    st.reason = strtrim(['recovered on pass 2; ' st.reason]);
                else
                    fprintf('  Pass 2 also failed: %s\n', st2.reason);
                    st.reason = [st.reason ' (pass 2 also failed)'];
                end
            end
            if st.ok
                nOK = nOK + 1; statusStr = 'OK'; reason = '';
                fprintf('\n');
            elseif st.partial
                nPartial = nPartial + 1; statusStr = 'PARTIAL'; reason = st.reason;
                fprintf('PARTIAL: %s\nReason: %s (tracked positions saved)\nContinuing...\n\n', it.name, reason);
            else
                nFailed = nFailed + 1; statusStr = 'FAILED'; reason = st.reason;
                fprintf('FAILED: %s\nReason: %s\nContinuing to next video...\n\n', it.name, reason);
            end
            r = row_from_item(it, opts.batchLabel, statusStr, reason, st.firstValidFrame, st.fps, resultsDir);
        catch ME
            nFailed = nFailed + 1;
            r = row_from_item(it, opts.batchLabel, 'FAILED', ME.message, NaN, NaN, resultsDir);
            fprintf('FAILED: %s\nReason: %s\nContinuing to next video...\n\n', it.name, ME.message);
        end
        log_row(L, i, total, r);
    end

    log_finalize(L, struct('total',total,'nOK',nOK,'nPartial',nPartial, ...
                           'nFailed',nFailed,'nSkipped',nSkipped));
    fprintf('\n=== BATCH COMPLETE ===\n');
    fprintf('OK=%d  PARTIAL=%d  FAILED=%d  SKIPPED=%d  (of %d)\n', ...
            nOK, nPartial, nFailed, nSkipped, total);
    fprintf('Log        : %s\n', L.txt);
    fprintf('Progress   : %s\n', L.csv);
    fprintf('Retry list : %s\n', L.retry);
end

% ═════════════════════════════════════════════════════════════════════════
%  MODE 3 — RERUN / DEBUG
% ═════════════════════════════════════════════════════════════════════════
function run_rerun(CFG, target, opts)
    clc; fprintf('\n=== RERUN SELECTED ===\n');
    if nargin < 3 || ~isstruct(opts), opts = normalize_opts(struct()); end
    if isempty(target)
        [file, inDir] = uigetfile({'*.avi;*.AVI;*.mp4;*.MP4;*.mov;*.MOV','Video Files'}, ...
            'Select the video for the trial to rerun');
        if isequal(file, 0), fprintf('Cancelled.\n'); return; end
        target = fullfile(inDir, file);
    end
    it = infer_item(target);
    if ~it.ok, error('%s', it.reason); end

    fields = inputdlg( ...
        {'Material (GB/CHIN)','Batch (blank=none)','Drop height (mm)', ...
         'Trial number','Container (full/shallow)', ...
         'Model (blank=none; Default/Tight/Wide Model)'}, ...
        'Confirm trial info', 1, ...
        {it.material, '', num2str(it.dropHeight_mm), num2str(it.trialNum), ...
         it.container, opts.model});
    if isempty(fields), fprintf('Cancelled.\n'); return; end
    m = build_meta(fields{1}, normalize_batch_label(fields{2}), ...
                   str2double(fields{3}), str2double(fields{4}), ...
                   fields{5}, fields{6});

    outRoot = resolve_output_root(CFG.outputRoot);
    if isempty(outRoot), fprintf('Cancelled.\n'); return; end
    [framesDir, detDir, resultsDir] = build_leaf_dirs(outRoot, m);

    stage = menu('Rerun stage', ...
        'Re-track only (reuse frames + detections)', ...
        'Re-detect + track (reuse frames)', ...
        'Full (re-export everything)');
    if stage == 0, fprintf('Cancelled.\n'); return; end

    cfg = base_cfg(CFG, m, target, framesDir, detDir, resultsDir);
    cfg.interactive = true; cfg.makeQA = true;
    switch stage
        case 1, cfg.reuseFrames = true;  cfg.reuseDetections = true;
        case 2, cfg.reuseFrames = true;  cfg.reuseDetections = false;
        case 3, cfg.reuseFrames = false; cfg.reuseDetections = false;
    end

    try
        st = process_one_trial(cfg);
        if ~st.ok && st.partial
            fprintf('PARTIAL: %s\nReason: %s (raw tracks saved)\n', it.name, st.reason);
        elseif ~st.ok
            fprintf('FAILED: %s\nReason: %s\n', it.name, st.reason);
        end
    catch ME
        fprintf('FAILED: %s\nReason: %s\n', it.name, ME.message);
    end
end

% ═════════════════════════════════════════════════════════════════════════
%  HELPERS
% ═════════════════════════════════════════════════════════════════════════
function label = normalize_batch_label(label)
%NORMALIZE_BATCH_LABEL  A bare number means "Batch <n>".
%   The batch label becomes a directory level, so '5' and 'Batch 5' would
%   otherwise create two sibling trees for the same batch and split its trials
%   across both. Anything that is not purely digits is left exactly as typed,
%   and empty stays empty.
    label = strtrim(char(string(label)));
    if isempty(label), label = ''; return; end
    if ~isempty(regexp(label, '^\d+$', 'once'))
        newLabel = ['Batch ' label];
        fprintf('Batch label normalized: "%s" -> "%s"\n', label, newLabel);
        label = newLabel;
    end
end

function opts = normalize_opts(opts)
    if isempty(opts) || ~isstruct(opts), opts = struct(); end
    if ~isfield(opts,'dryRun')  || isempty(opts.dryRun),  opts.dryRun  = false; end
    if ~isfield(opts,'limit')   || isempty(opts.limit),   opts.limit   = 0;     end
    if ~isfield(opts,'batchLabel'),                       opts.batchLabel = ''; end
    opts.batchLabel = normalize_batch_label(opts.batchLabel);
    if ~isfield(opts,'policy')  || isempty(opts.policy),  opts.policy  = 'reuse'; end
    if ~isfield(opts,'backupParams')
        % Pass-2 recovery: only used when pass 1 finds no valid frame.
        opts.backupParams = struct('sensitivity',0.92, 'edgeThresh',0.06, ...
                                   'radiusRange',[7 36]);
    end
    if ~isfield(opts,'videoListFile'),                    opts.videoListFile = ''; end
    if ~isfield(opts,'model')   || isempty(opts.model),   opts.model   = ''; end
    % autoWindow narrows the EXPORT/DETECT range only. It never changes which
    % trials are selected, and it is independent of dryRun and policy.
    if ~isfield(opts,'autoWindow') || isempty(opts.autoWindow)
        opts.autoWindow = true;
    end
    if ~isfield(opts,'windowPad') || isempty(opts.windowPad)
        opts.windowPad = [200 500];
    end
    opts.windowPad = double(opts.windowPad(:).');
    if numel(opts.windowPad) == 1, opts.windowPad = opts.windowPad([1 1]); end
    opts.policy = lower(char(opts.policy));
    opts.model  = strtrim(char(opts.model));
end

function cfg = apply_auto_window(cfg, opts, tag)
%APPLY_AUTO_WINDOW  Attach the export/detect window to cfg, if enabled.
%   Narrows the EXPORT/DETECT range only. Trial selection, dryRun and policy
%   are untouched. On any pre-scan failure the window is left unset, which
%   process_one_trial reads as the full range.
    cfg.autoWindow = isequal(opts.autoWindow, true);
    if ~cfg.autoWindow, return; end
    try
        [ws, we, info] = auto_window(cfg.videoPath, opts.windowPad, tag);
        cfg.windowStart = ws;
        cfg.windowEnd   = we;
        fprintf('autoWindow: frames %d-%d of %d (%.1f%%)\n', ...
                ws, we, info.nFrames, info.fracKept);
    catch ME
        cfg.autoWindow  = false;
        cfg.windowStart = [];
        cfg.windowEnd   = [];
        warning('process_trial:autoWindowFailed', ...
            ['%s: pre-scan failed (%s). Falling back to the FULL frame range ' ...
             'so nothing is silently skipped.'], tag, ME.message);
    end
end

% ═════════════════════════════════════════════════════════════════════════
function [winStart, winEnd, info] = auto_window(videoPath, pad, tag)
%AUTO_WINDOW  Frame range worth exporting, from a cheap red-presence pre-scan.
%
%   New-campaign clips run ~12 s (~40k frames) around a ~2k-frame event, so
%   exporting and detecting every frame spends almost all of its time on empty
%   bed. This walks the video once with VideoReader -- no PNG export, no
%   imfindcircles -- and flags each frame that contains any red pixel at all:
%
%       any(R > 150 & G < 100, 'all')
%
%   the same crude mask diag_raw_clip and survey_capture_integrity use. It is
%   deliberately not a detector: it only has to bracket the markers.
%
%   Window = [firstRed - pad(1), lastRed + pad(2)], clamped to [1, nFrames].
%   pad(1) buys pre-impact context for the velocity-peak search; pad(2) buys
%   the settling tail find_stop extrapolates across.
%
%   If NO frame is flagged the full range is returned with a NO_RED_CONTENT
%   warning, so a clip is never silently narrowed to nothing.

    winStart = NaN; winEnd = NaN;
    info = struct('ok',false,'reason','','nFrames',NaN, ...
                  'firstRed',NaN,'lastRed',NaN,'fracKept',NaN);

    v  = VideoReader(videoPath);
    nF = floor(v.Duration * v.FrameRate);
    info.nFrames = nF;

    firstRed = NaN; lastRed = NaN; k = 0;
    while hasFrame(v)
        k = k + 1;
        f = readFrame(v);
        if any(f(:,:,1) > 150 & f(:,:,2) < 100, 'all')
            if isnan(firstRed), firstRed = k; end
            lastRed = k;
        end
    end
    if k > 0, nF = max(nF, k); info.nFrames = nF; end

    if isnan(firstRed)
        winStart = 1; winEnd = nF;
        info.reason   = 'NO_RED_CONTENT';
        info.fracKept = 100;
        warning('process_trial:noRedContent', ...
            ['%s: no frame contains red marker pixels, so the impact window ' ...
             'could not be located. Falling back to the FULL range (1-%d) so ' ...
             'nothing is silently skipped. Check lighting, focus, or the clip ' ...
             'itself with scripts/diag_raw_clip.m.'], tag, nF);
        return
    end

    info.firstRed = firstRed;
    info.lastRed  = lastRed;
    winStart = max(1,  firstRed - pad(1));
    winEnd   = min(nF, lastRed  + pad(2));
    info.ok  = true;
    info.fracKept = 100 * (winEnd - winStart + 1) / max(nF,1);
end

function md = normalize_mode(mode)
    if isnumeric(mode), map = {'single','batch','rerun'}; md = map{mode}; return; end
    switch lower(strtrim(char(mode)))
        case {'1','single','manual'}, md = 'single';
        case {'2','batch'},           md = 'batch';
        case {'3','rerun','debug'},   md = 'rerun';
        otherwise,                    md = lower(strtrim(char(mode)));
    end
end

function params = default_detect_params()
    params = struct();
    params.radiusRange       = [9 33];
    params.sensitivity       = 0.85;
    params.edgeThresh        = 0.10;
    params.polarity          = 'bright';
    params.alphaG            = 0.50;
    params.betaB             = 0.50;
    params.doCLAHE           = true;
    params.medianK           = 3;
    params.showPreviewEveryN = 100;
    params.showCenters       = true;
    params.heightLabel       = '';
end

function m = build_meta(material, batchName, dropHeight_mm, trialNum, container, model)
    % model is optional. When '' the tag is exactly what it was before models
    % existed, so previously processed trials keep resolving to the same names.
    if nargin < 6 || isempty(model), model = ''; end
    container = lower(strtrim(container));
    if isempty(container), container = 'full'; end
    m = struct();
    m.material        = upper(strtrim(material));
    m.batchName       = strtrim(batchName);
    m.dropHeight_mm   = dropHeight_mm;
    m.trialNum        = trialNum;
    m.container       = container;
    m.model           = strtrim(char(model));
    m.heightLabel     = sprintf('%dmm', round(dropHeight_mm));
    m.trialParent     = sprintf('%dmm_T%02d', round(dropHeight_mm), trialNum);
    if isempty(m.model)
        m.trialTag = sprintf('%s_%s', m.trialParent, container);
    else
        % 'Tight Model' -> 'tight'. Appended so the same drop height and trial
        % number in a different model cannot collide on disk.
        suffix     = lower(char(erase(string(m.model), [" Model","Model"," ","_"])));
        m.trialTag = sprintf('%s_%s_%s', m.trialParent, container, suffix);
    end
    m.firstValidFrame = NaN;
    m.fps_true        = NaN;
    % Absolute video-frame span of the exported set; filled in by
    % process_one_trial once the window is known. Every frame index downstream
    % is relative to windowStart.
    m.windowStart     = NaN;
    m.windowEnd       = NaN;
    m.autoWindow      = false;
end

function cfg = base_cfg(CFG, m, videoPath, framesDir, detDir, resultsDir)
    cfg = struct();
    cfg.videoPath        = videoPath;
    cfg.framesDir        = framesDir;
    cfg.detDir           = detDir;
    cfg.resultsDir       = resultsDir;
    cfg.meta             = m;
    cfg.params           = CFG.params;
    % Per-model calibration: bed points and impactDistPx were measured per foot.
    % Falls back to the single global calibration when no model is set.
    if isfield(m,'model') && ~isempty(m.model)
        cfg.calib = get_calibration_model(m.model, m.dropHeight_mm, m.container);
    else
        cfg.calib = CFG.calib;
    end
    cfg.nExpectedMarkers = CFG.nExpectedMarkers;
    cfg.filterType       = CFG.filterType;
    cfg.interactive      = false;
    cfg.makeQA           = false;
    cfg.reuseFrames      = false;
    cfg.reuseDetections  = false;
end

function [framesDir, detDir, resultsDir] = build_leaf_dirs(outputRoot, m)
    % Model level sits directly under the batch level, so a tree written with
    % no model is byte-identical to the pre-model layout:
    %   with model : GB\Batch 5\Tight Model\165mm_T03\full
    %   without    : GB\Batch 5\165mm_T03\full
    hasModel = isfield(m,'model') && ~isempty(m.model);
    if isempty(m.batchName)
        if hasModel
            rel = fullfile(m.material, m.model, m.trialParent, m.container);
        else
            rel = fullfile(m.material, m.trialParent, m.container);
        end
    else
        if hasModel
            rel = fullfile(m.material, m.batchName, m.model, m.trialParent, m.container);
        else
            rel = fullfile(m.material, m.batchName, m.trialParent, m.container);
        end
    end
    % Exported frames are a DISPOSABLE intermediate (regenerable from the .avi).
    % To keep millions of PNGs out of Dropbox, they go under a LOCAL scratch root
    % if one is set via the JERBOA_FRAMES_ROOT environment variable; otherwise
    % they fall back to outputRoot (original behaviour). Detections and results
    % — the small, precious outputs — always stay under outputRoot.
    framesRoot = getenv('JERBOA_FRAMES_ROOT');
    if isempty(framesRoot), framesRoot = outputRoot; end
    framesDir  = fullfile(framesRoot, '01_FRAMES',           rel);
    detDir     = fullfile(outputRoot, '02_SAVED_DETECTIONS', rel);
    resultsDir = fullfile(outputRoot, '03_RESULTS',          rel);
end

% resolve_output_root moved to src/resolve_output_root.m (shared with track_tracers_2)

function items = items_from_list(listFile)
    txt   = fileread(listFile);
    lines = regexp(strtrim(txt), '\r?\n', 'split');
    lines = lines(~cellfun('isempty', lines));
    items = repmat(infer_item('seed'), 0, 1);
    for k = 1:numel(lines)
        items(k,1) = infer_item(strtrim(lines{k})); %#ok<AGROW>
    end
end

% ── row helpers ───────────────────────────────────────────────────────────
function r = row_from_item(it, batchLabel, status, reason, fvf, fps, resultsDir)
    r = struct('name',it.name, 'fullpath',it.fullpath, 'material',it.material, ...
        'batch',batchLabel, 'dropHeight_mm',it.dropHeight_mm, 'trialNum',it.trialNum, ...
        'container',it.container, 'firstValidFrame',fvf, 'fps',fps, ...
        'status',status, 'reason',reason, 'resultsDir',resultsDir);
end

% ── dry-run report (adapter → shared src/write_dryrun_report.m) ───────────
function print_and_save_dryrun(items, outRoot, batchLabel, inputDesc, model)
    if nargin < 5, model = ''; end
    rows = repmat(struct('head','','ok',false,'pathLines',{{}}), numel(items), 1);
    for i = 1:numel(items)
        it = items(i);
        if it.ok
            m = build_meta(it.material, batchLabel, it.dropHeight_mm, it.trialNum, ...
                           it.container, model);
            [fr, de, re] = build_leaf_dirs(outRoot, m);
            rows(i).head = sprintf(['[%3d] OK   %-26s  mat=%s cont=%s drop=%gmm ' ...
                                    'trial=T%02d model=%s tag=%s'], ...
                i, it.name, it.material, it.container, it.dropHeight_mm, ...
                it.trialNum, local_dash(m.model), m.trialTag);
            rows(i).ok        = true;
            rows(i).pathLines = {'frames', fr; 'detect', de; 'results', re};
        else
            rows(i).head = sprintf('[%3d] SKIP %-26s  reason: %s', i, it.name, it.reason);
            rows(i).ok   = false;
        end
    end
    write_dryrun_report(outRoot, inputDesc, rows);
end

% ── crash-safe incremental logger (adapters → shared src/batch_log_*.m) ────
function L = log_init(logDir, inputDesc, total, opts)
    stamp  = datestr(now, 'yyyymmdd_HHMMSS');
    header = sprintf([ ...
        '============================================================\n' ...
        ' BATCH LOG  %s\n' ...
        ' Input  : %s\n' ...
        ' Policy : %s   Limit : %d   Batch : "%s"\n' ...
        ' Model  : %s\n' ...
        ' Total videos : %d\n' ...
        '============================================================\n\n'], ...
        stamp, inputDesc, opts.policy, opts.limit, opts.batchLabel, ...
        local_dash(opts.model), total);
    csvHeader = sprintf(['idx,name,material,batch,dropHeight_mm,trialNum,container,', ...
                         'status,firstValidFrame,fps,reason,resultsDir,fullpath\n']);
    L = batch_log_init(logDir, stamp, ...
        {'batch_log_','batch_progress_','retry_failed_'}, header, csvHeader);
end

function log_row(L, i, total, r)
    txt = sprintf('[%3d/%3d] %-26s  %s\n', i, total, r.name, r.status);
    txt = [txt sprintf('          material=%s  batch=%s  drop=%gmm  trial=%s  container=%s\n', ...
        r.material, r.batch, r.dropHeight_mm, ttag(r.trialNum), r.container)];
    txt = [txt sprintf('          firstValidFrame=%s  fps=%s\n', ...
        num2str(r.firstValidFrame), fmtn(r.fps))];
    if ~isempty(r.reason)
        txt = [txt sprintf('          note   : %s\n', r.reason)];
    end
    if any(strcmpi(r.status,{'OK','PARTIAL'}))
        txt = [txt sprintf('          results: %s\n', r.resultsDir)];
    end
    txt = [txt sprintf('\n')];

    csvRow = sprintf('%d,%s,%s,%s,%g,%s,%s,%s,%s,%s,"%s",%s,%s\n', ...
        i, r.name, r.material, r.batch, r.dropHeight_mm, ...
        num2str(r.trialNum), r.container, r.status, ...
        num2str(r.firstValidFrame), fmtn(r.fps), r.reason, r.resultsDir, r.fullpath);

    retryLine = '';
    if any(strcmpi(r.status, {'FAILED','PARTIAL'})) && ~isempty(r.fullpath)
        retryLine = r.fullpath;
    end
    batch_log_row(L, txt, csvRow, retryLine);
end

function log_finalize(L, summ)
    txt = sprintf([ ...
        '============================================================\n' ...
        ' DONE  total=%d  OK=%d  PARTIAL=%d  FAILED=%d  SKIPPED=%d\n' ...
        '============================================================\n'], ...
        summ.total, summ.nOK, summ.nPartial, summ.nFailed, summ.nSkipped);
    batch_log_finalize(L, txt);
end

function s = ttag(n)
    if isnan(n), s = 'T??'; else, s = sprintf('T%02d', n); end
end
function s = fmtn(x)
    if isnan(x), s = 'NaN'; else, s = sprintf('%.4f', x); end
end
function s = local_dash(x)
    if isempty(x), s = '(none)'; else, s = char(x); end
end

function [material, batchName] = infer_material_batch(someDir) %#ok<DEFNU>
    d = char(someDir);
    if ~isempty(d) && (d(end)==filesep || d(end)=='/'), d(end) = []; end
    [parent, batchName] = fileparts(d);
    [~, material]       = fileparts(parent);
    material = upper(material);
end