function T = preflight_new_models(varargin)
%PREFLIGHT_NEW_MODELS  Confirm the calibration and output path for every new
%   trial BEFORE anything is written. READ-ONLY.
%
%   Builds the exact meta, tag, calibration and results path each video would
%   receive, then reports collisions against what is already on disk. Nothing
%   is processed and nothing is written.
%
%   USAGE
%       T = preflight_new_models;                    % all three models
%       T = preflight_new_models('ZeroDropOnly',true);
%       writetable(T, 'F:\preflight.csv');           % if you want a record
%
%   OPTIONS
%       'RawRoot'       default F:\ME_GRANULAB\Test Batches\Batch 5
%       'Root'          default F:\ME_GRANULAB\JerboaImpact
%       'BatchLabel'    default 'Batch 5'
%       'Models'        default all three
%       'ZeroDropOnly'  only h = 0 trials (default false)

opt.RawRoot      = 'F:\ME_GRANULAB\Test Batches\Batch 5';
opt.Root         = 'F:\ME_GRANULAB\JerboaImpact';
opt.BatchLabel   = 'Batch 5';
opt.Models       = ["Default Model","Tight Model","Wide Model"];
opt.ZeroDropOnly = false;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

fprintf('\n=== preflight_new_models  (READ-ONLY) ===\n');

E = dir(fullfile(opt.Root,'03_RESULTS','**','*_tracks.mat'));
existing = string(erase({E.name},'_tracks.mat'))';

rows = {};
for m = 1:numel(opt.Models)
    model = opt.Models(m);
    fdir  = fullfile(opt.RawRoot, char(model), 'GB', 'full');
    if ~isfolder(fdir), fprintf('MISSING: %s\n', fdir); continue; end
    V = dir(fullfile(fdir,'*.avi')); V = V(~[V.isdir]);
    for i = 1:numel(V)
        [~, stem] = fileparts(V(i).name);
        tok = regexp(stem,'^(\d+)mm_T(\d+)$','tokens','once');
        if isempty(tok), continue; end
        h = str2double(tok{1});  tn = str2double(tok{2});
        if opt.ZeroDropOnly && h ~= 0, continue; end

        % model suffix on the tag, model level in the path
        suffix   = lower(erase(model," Model"));
        trialTag = sprintf('%dmm_T%02d_full_%s', h, tn, suffix);
        resDir   = fullfile(opt.Root,'03_RESULTS','GB',opt.BatchLabel, ...
                            char(model), sprintf('%dmm_T%02d',h,tn), 'full');

        c = get_calibration_model(model, h, 'full');
        rows{end+1} = table(model, string(V(i).name), h, tn, string(trialTag), ...
            c.impactDistPx, c.bedPoint2(2), h==0, ...
            ismember(string(trialTag), existing), string(resDir), ...
            'VariableNames',{'model','file','dropHeight_mm','trialNum','trialTag', ...
                'impactDistPx','bedPoint2_y','isZeroDrop','tagExists','resultsDir'}); %#ok<AGROW>
    end
end
if isempty(rows), error('preflight_new_models:none','No parseable videos found.'); end
T = vertcat(rows{:});

fprintf('\n--- calibration actually assigned, by group ---\n');
[g, gm, gh] = findgroups(T.model, T.bedPoint2_y);
for k = 1:max(g)
    s = T(g==k,:);
    fprintf('  %-14s bed (4,0)-(4,%2d)  trigger %9.3f  n=%3d  heights: %s\n', ...
        gm(k), gh(k), s.impactDistPx(1), height(s), ...
        strjoin(compose('%g', unique(s.dropHeight_mm)'), ', '));
end

fprintf('\n--- counts ---\n');
for mm = unique(T.model,'stable')'
    s = T(T.model==mm,:);
    fprintf('  %-14s %3d total | %2d zero-drop | %3d new tags | %3d tags already on disk\n', ...
        mm, height(s), sum(s.isZeroDrop), sum(~s.tagExists), sum(s.tagExists));
end

if any(T.tagExists)
    fprintf(['\n  %d tag(s) already exist. With the model suffix these should be 0 --\n' ...
             '  if not, the suffix is not disambiguating and processing would overwrite.\n'], ...
             sum(T.tagExists));
    disp(T(T.tagExists, {'model','file','trialTag'}));
else
    fprintf('\n  No tag collisions: every new trial writes to a fresh path.\n');
end

fprintf('\nNothing written. Review the table, then process.\n\n');
end