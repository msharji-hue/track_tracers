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
%     .policy        'reuse'(default) | 'resume' | 'overwrite'
%                      reuse     : reuse frames+detections if present, redo tracks
%                      resume    : additionally SKIP trials already having tracks
%                      overwrite : re-export, re-detect, re-track everything
%     .videoListFile a .txt of full video paths to process (e.g. a retry list)

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
        case 'single', run_single(CFG);
        case 'batch',  run_batch(CFG, inputTarget, opts);
        case 'rerun',  run_rerun(CFG, inputTarget);
        otherwise,     error('Unknown mode "%s" (single|batch|rerun)', mode);
    end
end

% ═════════════════════════════════════════════════════════════════════════
%  MODE 1 — SINGLE / MANUAL
% ═════════════════════════════════════════════════════════════════════════
function run_single(CFG)
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
         'Trial number','Container (full/shallow)'}, ...
        'Confirm trial info', 1, {it.material, '', defH, defT, it.container});
    if isempty(fields), fprintf('Cancelled.\n'); return; end

    m = build_meta(fields{1}, fields{2}, str2double(fields{3}), str2double(fields{4}), fields{5});
    outRoot = resolve_output_root(CFG.outputRoot);
    if isempty(outRoot), fprintf('Cancelled.\n'); return; end
    [framesDir, detDir, resultsDir] = build_leaf_dirs(outRoot, m);

    cfg = base_cfg(CFG, m, fullfile(inDir, file), framesDir, detDir, resultsDir);
    cfg.interactive = true; cfg.makeFigures = true;

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
        a = inputdlg({'Batch label (blank = none in path)','Limit (0 = all)'}, ...
                     'Batch options', 1, {opts.batchLabel, num2str(opts.limit)});
        if isempty(a), fprintf('Cancelled.\n'); return; end
        opts.batchLabel = strtrim(a{1});
        opts.limit      = max(0, round(str2double(a{2})));

        pol = menu('Redo policy', ...
            'reuse  : keep frames+detections, recompute tracks (default)', ...
            'resume : skip trials already done', ...
            'overwrite : redo everything');
        if pol == 0, fprintf('Cancelled.\n'); return; end
        pols = {'reuse','resume','overwrite'}; opts.policy = pols{pol};

        dr = menu('Run type', 'Dry run (list only, no processing)', 'Process for real');
        if dr == 0, fprintf('Cancelled.\n'); return; end
        opts.dryRun = (dr == 1);
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
        print_and_save_dryrun(items, outRoot, opts.batchLabel, inputDesc);
        return;
    end

    % ── Process (crash-safe incremental log) ─────────────────────────────
    total = numel(items);
    fprintf('Processing %d videos  (policy=%s, batch="%s")\n\n', total, opts.policy, opts.batchLabel);
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

        m = build_meta(it.material, opts.batchLabel, it.dropHeight_mm, it.trialNum, it.container);
        [framesDir, detDir, resultsDir] = build_leaf_dirs(outRoot, m);

        if strcmp(opts.policy,'resume') && ...
           exist(fullfile(resultsDir,'tracks',[m.trialTag '_tracks.mat']),'file')
            nSkipped = nSkipped + 1;
            r = row_from_item(it, opts.batchLabel, 'SKIPPED', 'already processed', NaN, NaN, resultsDir);
            fprintf('SKIPPED (already done): %s\n\n', it.name);
            log_row(L, i, total, r); continue;
        end

        cfg = base_cfg(CFG, m, it.fullpath, framesDir, detDir, resultsDir);
        cfg.interactive = false; cfg.makeFigures = true;
        if strcmp(opts.policy,'overwrite')
            cfg.reuseFrames = false; cfg.reuseDetections = false;
        else
            cfg.reuseFrames = true;  cfg.reuseDetections = true;
        end

        try
            st = process_one_trial(cfg);
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
function run_rerun(CFG, target)
    clc; fprintf('\n=== RERUN SELECTED ===\n');
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
         'Trial number','Container (full/shallow)'}, ...
        'Confirm trial info', 1, ...
        {it.material, '', num2str(it.dropHeight_mm), num2str(it.trialNum), it.container});
    if isempty(fields), fprintf('Cancelled.\n'); return; end
    m = build_meta(fields{1}, fields{2}, str2double(fields{3}), str2double(fields{4}), fields{5});

    outRoot = resolve_output_root(CFG.outputRoot);
    if isempty(outRoot), fprintf('Cancelled.\n'); return; end
    [framesDir, detDir, resultsDir] = build_leaf_dirs(outRoot, m);

    stage = menu('Rerun stage', ...
        'Re-track only (reuse frames + detections)', ...
        'Re-detect + track (reuse frames)', ...
        'Full (re-export everything)');
    if stage == 0, fprintf('Cancelled.\n'); return; end

    cfg = base_cfg(CFG, m, target, framesDir, detDir, resultsDir);
    cfg.interactive = true; cfg.makeFigures = true;
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
function opts = normalize_opts(opts)
    if isempty(opts) || ~isstruct(opts), opts = struct(); end
    if ~isfield(opts,'dryRun')  || isempty(opts.dryRun),  opts.dryRun  = false; end
    if ~isfield(opts,'limit')   || isempty(opts.limit),   opts.limit   = 0;     end
    if ~isfield(opts,'batchLabel'),                       opts.batchLabel = ''; end
    if ~isfield(opts,'policy')  || isempty(opts.policy),  opts.policy  = 'reuse'; end
    if ~isfield(opts,'videoListFile'),                    opts.videoListFile = ''; end
    opts.policy = lower(char(opts.policy));
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

function m = build_meta(material, batchName, dropHeight_mm, trialNum, container)
    container = lower(strtrim(container));
    if isempty(container), container = 'full'; end
    m = struct();
    m.material        = upper(strtrim(material));
    m.batchName       = strtrim(batchName);
    m.dropHeight_mm   = dropHeight_mm;
    m.trialNum        = trialNum;
    m.container       = container;
    m.heightLabel     = sprintf('%dmm', round(dropHeight_mm));
    m.trialParent     = sprintf('%dmm_T%02d', round(dropHeight_mm), trialNum);
    m.trialTag        = sprintf('%s_%s', m.trialParent, container);
    m.firstValidFrame = NaN;
    m.fps_true        = NaN;
end

function cfg = base_cfg(CFG, m, videoPath, framesDir, detDir, resultsDir)
    cfg = struct();
    cfg.videoPath        = videoPath;
    cfg.framesDir        = framesDir;
    cfg.detDir           = detDir;
    cfg.resultsDir       = resultsDir;
    cfg.meta             = m;
    cfg.params           = CFG.params;
    cfg.calib            = CFG.calib;
    cfg.nExpectedMarkers = CFG.nExpectedMarkers;
    cfg.filterType       = CFG.filterType;
    cfg.interactive      = false;
    cfg.makeFigures      = false;
    cfg.reuseFrames      = false;
    cfg.reuseDetections  = false;
end

function [framesDir, detDir, resultsDir] = build_leaf_dirs(outputRoot, m)
    if isempty(m.batchName)
        rel = fullfile(m.material, m.trialParent, m.container);
    else
        rel = fullfile(m.material, m.batchName, m.trialParent, m.container);
    end
    framesDir  = fullfile(outputRoot, '01_FRAMES',           rel);
    detDir     = fullfile(outputRoot, '02_SAVED_DETECTIONS', rel);
    resultsDir = fullfile(outputRoot, '03_RESULTS',          rel);
end

function outRoot = resolve_output_root(preset)
    if ~isempty(preset), outRoot = preset; return; end
    outRoot = uigetdir(pwd, ...
        'Select OUTPUT root (parent of 01_FRAMES / 02_SAVED_DETECTIONS / 03_RESULTS)');
    if isequal(outRoot, 0), outRoot = ''; end
end

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

% ── dry-run report ────────────────────────────────────────────────────────
function print_and_save_dryrun(items, outRoot, batchLabel, inputDesc)
    stamp  = datestr(now, 'yyyymmdd_HHMMSS');
    repDir = fullfile(outRoot, '03_RESULTS', '_batch_logs');
    if ~exist(repDir, 'dir'), mkdir(repDir); end
    repPath = fullfile(repDir, sprintf('dryrun_report_%s.txt', stamp));
    fid = fopen(repPath, 'w');

    fprintf('DRY RUN — %d videos found under %s\n\n', numel(items), inputDesc);
    fprintf(fid, 'DRY RUN  %s\nInput: %s\nVideos found: %d\n\n', stamp, inputDesc, numel(items));

    nOK = 0; nBad = 0;
    for i = 1:numel(items)
        it = items(i);
        if it.ok
            nOK = nOK + 1;
            m = build_meta(it.material, batchLabel, it.dropHeight_mm, it.trialNum, it.container);
            [fr, de, re] = build_leaf_dirs(outRoot, m);
            head = sprintf('[%3d] OK   %-26s  mat=%s cont=%s drop=%gmm trial=T%02d', ...
                i, it.name, it.material, it.container, it.dropHeight_mm, it.trialNum);
            fprintf('%s\n', head);
            fprintf(fid, '%s\n      frames : %s\n      detect : %s\n      results: %s\n\n', head, fr, de, re);
        else
            nBad = nBad + 1;
            head = sprintf('[%3d] SKIP %-26s  reason: %s', i, it.name, it.reason);
            fprintf('%s\n', head);
            fprintf(fid, '%s\n\n', head);
        end
    end
    fprintf('\nDry run: %d processable, %d unparseable.\nReport: %s\n', nOK, nBad, repPath);
    fprintf(fid, 'SUMMARY: %d processable, %d unparseable\n', nOK, nBad);
    fclose(fid);
end

% ── crash-safe incremental logger (writes per trial) ──────────────────────
function L = log_init(logDir, inputDesc, total, opts)
    if ~exist(logDir, 'dir'), mkdir(logDir); end
    L.stamp = datestr(now, 'yyyymmdd_HHMMSS');
    L.txt   = fullfile(logDir, sprintf('batch_log_%s.txt',      L.stamp));
    L.csv   = fullfile(logDir, sprintf('batch_progress_%s.csv', L.stamp));
    L.retry = fullfile(logDir, sprintf('retry_failed_%s.txt',   L.stamp));

    fid = fopen(L.txt, 'w');
    fprintf(fid, '============================================================\n');
    fprintf(fid, ' BATCH LOG  %s\n', L.stamp);
    fprintf(fid, ' Input  : %s\n', inputDesc);
    fprintf(fid, ' Policy : %s   Limit : %d   Batch : "%s"\n', opts.policy, opts.limit, opts.batchLabel);
    fprintf(fid, ' Total videos : %d\n', total);
    fprintf(fid, '============================================================\n\n');
    fclose(fid);

    fid = fopen(L.csv, 'w');
    fprintf(fid, ['idx,name,material,batch,dropHeight_mm,trialNum,container,', ...
                  'status,firstValidFrame,fps,reason,resultsDir,fullpath\n']);
    fclose(fid);
end

function log_row(L, i, total, r)
    fid = fopen(L.txt, 'a');     % append+close each trial => crash-safe
    fprintf(fid, '[%3d/%3d] %-26s  %s\n', i, total, r.name, r.status);
    fprintf(fid, '          material=%s  batch=%s  drop=%gmm  trial=%s  container=%s\n', ...
        r.material, r.batch, r.dropHeight_mm, ttag(r.trialNum), r.container);
    fprintf(fid, '          firstValidFrame=%s  fps=%s\n', ...
        num2str(r.firstValidFrame), fmtn(r.fps));
    if ~isempty(r.reason),     fprintf(fid, '          note   : %s\n', r.reason); end
    if any(strcmpi(r.status,{'OK','PARTIAL'}))
        fprintf(fid, '          results: %s\n', r.resultsDir);
    end
    fprintf(fid, '\n');
    fclose(fid);

    fid = fopen(L.csv, 'a');
    fprintf(fid, '%d,%s,%s,%s,%g,%s,%s,%s,%s,%s,"%s",%s,%s\n', ...
        i, r.name, r.material, r.batch, r.dropHeight_mm, ...
        num2str(r.trialNum), r.container, r.status, ...
        num2str(r.firstValidFrame), fmtn(r.fps), r.reason, r.resultsDir, r.fullpath);
    fclose(fid);

    if any(strcmpi(r.status, {'FAILED','PARTIAL'})) && ~isempty(r.fullpath)
        fid = fopen(L.retry, 'a'); fprintf(fid, '%s\n', r.fullpath); fclose(fid);
    end
end

function log_finalize(L, summ)
    fid = fopen(L.txt, 'a');
    fprintf(fid, '============================================================\n');
    fprintf(fid, ' DONE  total=%d  OK=%d  PARTIAL=%d  FAILED=%d  SKIPPED=%d\n', ...
        summ.total, summ.nOK, summ.nPartial, summ.nFailed, summ.nSkipped);
    fprintf(fid, '============================================================\n');
    fclose(fid);
end

function s = ttag(n)
    if isnan(n), s = 'T??'; else, s = sprintf('T%02d', n); end
end
function s = fmtn(x)
    if isnan(x), s = 'NaN'; else, s = sprintf('%.4f', x); end
end

function [material, batchName] = infer_material_batch(someDir) %#ok<DEFNU>
    d = char(someDir);
    if ~isempty(d) && (d(end)==filesep || d(end)=='/'), d(end) = []; end
    [parent, batchName] = fileparts(d);
    [~, material]       = fileparts(parent);
    material = upper(material);
end
