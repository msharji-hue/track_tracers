function run_new_models(varargin)
%RUN_NEW_MODELS  Process only the new Batch 5 GB/full trials, per model.
%
%   Wraps process_trial('batch', ...) once per model, passing the model
%   through opts so build_meta tags it, build_leaf_dirs paths it, and
%   base_cfg picks up get_calibration_model. Nothing else in the pipeline
%   changes.
%
%   REQUIRES the four edits in process_trial_edits.txt.
%
%   USAGE
%       run_new_models('DryRun',true)                  % list only, no writes
%       run_new_models('Models',"Default Model", 'ZeroDropOnly',true)
%       run_new_models                                 % process the default set
%
%   OPTIONS
%       'RawRoot'       default F:\ME_GRANULAB\Test Batches\Batch 5
%       'Root'          default F:\ME_GRANULAB\JerboaImpact
%       'BatchLabel'    default 'Batch 5'
%       'Models'        default all three
%       'ZeroDropOnly'  only h = 0 trials (default false)
%       'DryRun'        list and plan, process nothing (default false)
%       'Policy'        'resume' (default) | 'reuse' | 'overwrite'
%                       resume SKIPS anything already having tracks, so a
%                       re-run after an interruption costs nothing.
%
%   DEFAULT SET. Default Model contributes its ZERO-DROP trials only: its
%   non-zero trials are already processed under the old un-suffixed tags, and
%   re-running them would leave two copies of the same data under two names.
%   Tight and Wide contribute everything.

opt.RawRoot      = 'F:\ME_GRANULAB\Test Batches\Batch 5';
opt.Root         = 'F:\ME_GRANULAB\JerboaImpact';
opt.BatchLabel   = 'Batch 5';
opt.Models       = ["Default Model","Tight Model","Wide Model"];
opt.ZeroDropOnly = false;
opt.DryRun       = false;
opt.Policy       = 'resume';
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

fprintf('\n=== run_new_models ===\n');
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

    % Default Model: zero-drop only, unless overridden
    zeroOnly = opt.ZeroDropOnly || strcmpi(erase(model," Model"),"Default");

    V = dir(fullfile(fdir,'*.avi')); V = V(~[V.isdir]);
    keep = strings(0);
    for i = 1:numel(V)
        [~, stem] = fileparts(V(i).name);
        tok = regexp(stem,'^(\d+)mm_T(\d+)$','tokens','once');
        if isempty(tok), continue; end
        if zeroOnly && str2double(tok{1}) ~= 0, continue; end
        keep(end+1) = string(fullfile(V(i).folder, V(i).name)); %#ok<AGROW>
    end

    fprintf('\n--- %s ---\n', model);
    fprintf('  %d of %d videos selected%s\n', numel(keep), numel(V), ...
        local_tern(zeroOnly, '  (zero-drop only)', ''));
    if isempty(keep), continue; end

    % An explicit video list is safer than a folder scan: process_trial reads
    % it verbatim, so nothing outside this list can be touched.
    listFile = fullfile(listDir, sprintf('videolist_%s_%s.txt', ...
                        matlab.lang.makeValidName(char(model)), stamp));
    fid = fopen(listFile,'w');
    fprintf(fid, '%s\n', keep);
    fclose(fid);
    fprintf('  list: %s\n', listFile);

    o = struct();
    o.videoListFile = listFile;
    o.batchLabel    = opt.BatchLabel;
    o.model         = char(model);      % -> build_meta -> tag, path, calib
    o.policy        = opt.Policy;
    o.dryRun        = opt.DryRun;
    o.limit         = 0;

    process_trial('batch', '', opt.Root, o);
end

fprintf('\n=== run_new_models done ===\n');
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
