function C = run_models(varargin)

opt.RawRoot      = 'D:\ME_GRANULAB\Test Batches\Batch 5';
opt.Root         = 'D:\ME_GRANULAB\JerboaImpact';
opt.BatchLabel   = 'Batch 5';
opt.Models       = ["Tight Model","Wide Model"];
opt.Heights      = [];
opt.ZeroDropOnly = false;
opt.Policy       = 'retry';
opt.BackupParams = [];
opt.DryRun       = false;
opt.CoverageOnly = false;
opt.Save         = false;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

fprintf('\n=== run_models ===\n');

% ── 1) Process ────────────────────────────────────────────────────────────
if ~opt.CoverageOnly
    fprintf('policy: %s   dryRun: %d\n', opt.Policy, opt.DryRun);
    listDir = fullfile(opt.Root,'03_RESULTS','_batch_logs');
    if ~isfolder(listDir), mkdir(listDir); end
    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));

    for m = 1:numel(opt.Models)
        model = opt.Models(m);
        fdir  = fullfile(opt.RawRoot, char(model), 'GB', 'full');
        if ~isfolder(fdir)
            fprintf('\nSKIP %s -- folder not found:\n  %s\n', model, fdir);
            continue
        end
        rdir     = fullfile(opt.Root,'03_RESULTS','GB',opt.BatchLabel,char(model));
        zeroOnly = opt.ZeroDropOnly || strcmpi(erase(model," Model"),"Default");
        suffix   = lower(char(erase(model, [" Model","Model"," ","_"])));

        % tags that already tracked -- these are left alone
        Tk   = dir(fullfile(rdir,'**','*_tracks.mat'));
        done = string(erase({Tk.name},'_tracks.mat'));

        V = dir(fullfile(fdir,'*.avi')); V = V(~[V.isdir]);
        keep = strings(0); nSkip = 0;
        for i = 1:numel(V)
            [~, stem] = fileparts(V(i).name);
            tok = regexp(stem,'^(\d+)mm_T(\d+)$','tokens','once');
            if isempty(tok), continue; end
            h = str2double(tok{1});
            if zeroOnly && h ~= 0, continue; end
            if ~zeroOnly && ~isempty(opt.Heights) && ~ismember(h, opt.Heights), continue; end
            tag = sprintf('%s_full_%s', stem, suffix);
            if ismember(string(tag), done), nSkip = nSkip + 1; continue; end
            keep(end+1) = string(fullfile(V(i).folder, V(i).name)); %#ok<AGROW>
        end

        fprintf('\n--- %s ---\n', model);
        fprintf('  %d videos present | %d already tracked (skipped) | %d to process\n', ...
                numel(V), nSkip, numel(keep));
        if isempty(keep), continue; end

        % An explicit list is safer than a folder scan: process_trial reads it
        % verbatim, so nothing outside the list can be touched.
        listFile = fullfile(listDir, sprintf('videolist_%s_%s.txt', ...
                            matlab.lang.makeValidName(char(model)), stamp));
        fid = fopen(listFile,'w');
        fprintf(fid, '%s\n', keep);
        fclose(fid);
        fprintf('  list: %s\n', listFile);

        o = struct('videoListFile',listFile, 'batchLabel',opt.BatchLabel, ...
                   'model',char(model), 'policy',opt.Policy, ...
                   'dryRun',opt.DryRun, 'limit',0);
        if ~isempty(opt.BackupParams), o.backupParams = opt.BackupParams; end
        process_trial('batch', '', opt.Root, o);
    end
end

% ── 2) Coverage ───────────────────────────────────────────────────────────
C = coverage_table(opt);
if isempty(C), fprintf('\nNo videos found.\n\n'); return; end

fprintf('\n=== COVERAGE ===\n');
for mm = unique(C.model,'stable')'
    s = C(C.model == mm, :);

    z = s(s.dropHeight_mm == 0, :);
    if ~isempty(z)
        fprintf('\n%s  --  d0 (h = 0)\n', mm);
        fprintf('   %d videos | %d tracked | %d untracked | %.0f%%\n', ...
            z.nVideos(1), z.nTracked(1), z.nUntracked(1), 100*z.successRate(1));
    end

    d = s(s.dropHeight_mm > 0, :);
    if isempty(d), continue; end
    fprintf('\n%s  --  drop-height trials\n', mm);
    fprintf('   %-10s %8s %8s %10s %8s\n', 'height', 'videos', 'tracked', 'untracked', 'rate');
    for i = 1:height(d)
        fprintf('   %7g mm %8d %8d %10d %7.0f%%\n', d.dropHeight_mm(i), ...
            d.nVideos(i), d.nTracked(i), d.nUntracked(i), 100*d.successRate(i));
    end
    fprintf('   %-10s %8d %8d %10d %7.0f%%\n', 'TOTAL', sum(d.nVideos), ...
        sum(d.nTracked), sum(d.nUntracked), 100*sum(d.nTracked)/sum(d.nVideos));
end
fprintf('\n');

if opt.Save
    st = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    p  = fullfile(opt.Root,'03_RESULTS','_batch_logs', ...
                  sprintf('coverage_models_%s.csv', st));
    writetable(C, p);
    fprintf('wrote: %s\n\n', p);
end
end

% ─────────────────────────────────────────────────────────────────────────
function C = coverage_table(opt)
%COVERAGE_TABLE  Videos present vs trials tracked, per model and drop height.
%   Counts videos on disk against tracks on disk, so it reflects the current
%   state regardless of when or how a trial was processed.
rows = {};
for m = 1:numel(opt.Models)
    model  = opt.Models(m);
    fdir   = fullfile(opt.RawRoot, char(model), 'GB', 'full');
    rdir   = fullfile(opt.Root,'03_RESULTS','GB',opt.BatchLabel,char(model));
    if ~isfolder(fdir), continue; end
    suffix = lower(char(erase(model, [" Model","Model"," ","_"])));

    Tk   = dir(fullfile(rdir,'**','*_tracks.mat'));
    done = string(erase({Tk.name},'_tracks.mat'));

    V = dir(fullfile(fdir,'*.avi')); V = V(~[V.isdir]);
    hv = nan(numel(V),1); tk = false(numel(V),1);
    for i = 1:numel(V)
        [~,stem] = fileparts(V(i).name);
        tok = regexp(stem,'^(\d+)mm_T(\d+)$','tokens','once');
        if isempty(tok), continue; end
        hv(i) = str2double(tok{1});
        tk(i) = ismember(string(sprintf('%s_full_%s', stem, suffix)), done);
    end

    for h = unique(hv(isfinite(hv)))'
        k = hv == h;
        rows{end+1} = table(model, h, sum(k), sum(k & tk), sum(k & ~tk), ...
            sum(k & tk)/sum(k), ...
            'VariableNames',{'model','dropHeight_mm','nVideos','nTracked', ...
                             'nUntracked','successRate'}); %#ok<AGROW>
    end
end
if isempty(rows), C = table(); else, C = vertcat(rows{:}); end
end
