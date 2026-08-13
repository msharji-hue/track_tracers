function T = review_model_trials(varargin)
%REVIEW_MODEL_TRIALS  One row per trial across the foot models. READ-ONLY.
%
%   Surveys every trial that produced detections, whether or not it tracked,
%   so failures and marginal successes appear side by side rather than only
%   the ones that happened to fail loudly.
%
%   The key columns are the detection-count summary. A trial that initialises
%   on a single transient frame is not a real success: nEightFrames = 1 means
%   the 8-marker condition was met exactly once, which is fragile.
%
%   USAGE
%       T = review_model_trials;
%       T = review_model_trials('Models', "Tight Model");
%       sortrows(T, 'nEightFrames')            % worst first
%       T(T.nEightFrames <= 2, :)              % trials on a knife edge
%
%   OPTIONS
%       'Root'    default F:\ME_GRANULAB\JerboaImpact
%       'Models'  default all three
%       'Save'    write a timestamped CSV to _batch_logs (default false)

opt.Root   = 'F:\ME_GRANULAB\JerboaImpact';
opt.Models = ["Default Model","Tight Model","Wide Model"];
opt.Save   = false;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

fprintf('\n=== review_model_trials  (READ-ONLY) ===\n');

rows = {};
for m = 1:numel(opt.Models)
    model  = opt.Models(m);
    detDir = fullfile(opt.Root,'02_SAVED_DETECTIONS','GB','Batch 5',char(model));
    resDir = fullfile(opt.Root,'03_RESULTS','GB','Batch 5',char(model));
    if ~isfolder(detDir), fprintf('  no detections for %s\n', model); continue; end

    D = dir(fullfile(detDir,'**','*_detections.mat'));
    fprintf('  %-14s %3d trials with detections\n', model, numel(D));

    for i = 1:numel(D)
        tag = string(erase(D(i).name,'_detections.mat'));
        S = load(fullfile(D(i).folder, D(i).name));
        n = S.det.detect.nDetected(:);

        % Did it track? Tracks exist only when an 8-marker frame was found.
        Tk = dir(fullfile(resDir,'**',char(tag) + "_tracks.mat"));
        tracked = ~isempty(Tk);
        fvf = NaN; nTr = NaN;
        if tracked
            K = load(fullfile(Tk(1).folder, Tk(1).name), 'tracks');
            fvf = K.tracks.firstValidFrame;
            nTr = K.tracks.nTracked;
        end

        % first index whose count is at the modal value, i.e. when the rod
        % has settled into view rather than sweeping through
        md = mode(n(n > 0));
        if isempty(md) || isnan(md)
            fEnter = NaN;
        else
            fEnter = find(n == md, 1);
        end

        h = regexp(char(tag), '^(\d+)mm', 'tokens', 'once');
        rows{end+1} = table(model, tag, str2double(h{1}), numel(n), ...
            sum(n == 0), md, sum(n == md), max(n), sum(n == 8), ...
            fEnter, tracked, fvf, nTr, ...
            'VariableNames',{'model','trialTag','dropHeight_mm','nFrames', ...
                'nZeroFrames','modalCount','nModalFrames','maxCount', ...
                'nEightFrames','firstNonZero','tracked','firstValidFrame', ...
                'nTracked'}); %#ok<AGROW>
    end
end
if isempty(rows), error('review_model_trials:none','No detections found.'); end
T = vertcat(rows{:});

% ------------------------------------------------------------- summary
fprintf('\n--- per model ---\n');
for mm = unique(T.model,'stable')'
    s = T(T.model == mm, :);
    fprintf('  %-14s %3d trials | %3d tracked (%.0f%%) | %3d never reached 8 markers\n', ...
        mm, height(s), sum(s.tracked), 100*mean(s.tracked), sum(s.nEightFrames == 0));
    fprintf('  %-14s modal marker count: %s\n', '', ...
        strjoin(compose('%d', unique(s.modalCount(~isnan(s.modalCount))))', ', '));
end

fprintf('\n--- how fragile are the successes? ---\n');
fprintf('  Frames meeting the 8-marker condition, among tracked trials:\n');
tr = T(T.tracked, :);
for k = [1 2 3 5 10]
    fprintf('    <= %2d frames : %3d of %3d\n', k, sum(tr.nEightFrames <= k), height(tr));
end
if any(tr.nEightFrames == 1)
    fprintf(['  %d tracked trial(s) met the condition on exactly ONE frame --\n' ...
             '  those initialised on a single transient and are not solid.\n'], ...
             sum(tr.nEightFrames == 1));
end

fprintf('\n--- trials worth looking at first (fewest 8-marker frames) ---\n');
disp(head(sortrows(T, 'nEightFrames'), 15));

if opt.Save
    st = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    p = fullfile(opt.Root,'03_RESULTS','_batch_logs', ...
                 sprintf('review_model_trials_%s.csv', st));
    writetable(T, p);
    fprintf('wrote: %s\n', p);
end
fprintf('\n');
end
