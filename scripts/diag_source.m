function diag_source(root, rawRoot, nPerModel)
% DIAG_SOURCE  Are the short Tight/Wide tracks a Stage-A truncation, or were
% the source videos themselves short?
%
%   For a sample of trials per model, compares the RAW video frame count to the
%   number of frames saved in *_tracks.mat.
%
%       rawFrames >> trackedFrames   -> Stage A truncated; reprocessing helps
%       rawFrames ~= trackedFrames   -> the capture was short; reprocessing
%                                       cannot recover the event
%
%   diag_source('D:\ME_GRANULAB\JerboaImpact', 'D:\ME_GRANULAB\Test Batches', 8)

    if nargin < 2, rawRoot = 'D:\ME_GRANULAB\Test Batches'; end
    if nargin < 3, nPerModel = 8; end

    % ---- index every raw video by basename ----
    exts = {'*.avi','*.mp4','*.mov','*.cine','*.mraw'};
    V = [];
    for e = 1:numel(exts)
        V = [V; dir(fullfile(rawRoot,'**',exts{e}))]; %#ok<AGROW>
    end
    fprintf('indexed %d raw video file(s) under %s\n', numel(V), rawRoot);
    if isempty(V)
        warning('No videos found - check rawRoot and the extension list.');
        return;
    end
    vname = strings(numel(V),1);
    for k = 1:numel(V), [~,b] = fileparts(V(k).name); vname(k) = lower(string(b)); end

    % ---- sample tracks per model ----
    D = dir(fullfile(root,'03_RESULTS','**','*_tracks.mat'));
    tg = strings(numel(D),1); mdl = tg;
    for k = 1:numel(D)
        t = string(erase(D(k).name,'_tracks.mat'));
        tg(k) = t;
        if endsWith(lower(t),"_tight"), mdl(k)="Tight";
        elseif endsWith(lower(t),"_wide"), mdl(k)="Wide";
        else, mdl(k)="Default"; end
    end

    rows = {};
    for M = ["Default","Tight","Wide"]
        idx = find(mdl==M);
        if isempty(idx), continue; end
        pick = idx(round(linspace(1, numel(idx), min(nPerModel,numel(idx)))));
        for j = pick(:)'
            tp = fullfile(D(j).folder, D(j).name);
            S = load(tp,'tracks','meta');
            nTracked = size(S.tracks.trackedX,2);
            fps = resolve_fps(tp, S.meta, S.tracks);

            % match a raw video: exact basename, else contains the tag
            hit = find(vname == lower(tg(j)), 1);
            if isempty(hit), hit = find(contains(vname, lower(tg(j))), 1); end
            if isempty(hit)
                rows(end+1,:) = {char(M), char(tg(j)), fps, nTracked, NaN, NaN, 'NO VIDEO MATCH'}; %#ok<AGROW>
                continue;
            end
            vp = fullfile(V(hit).folder, V(hit).name);
            try
                vr = VideoReader(vp);
                nRaw = vr.NumFrames;
            catch
                nRaw = NaN;
            end
            ratio = nRaw / nTracked;
            rows(end+1,:) = {char(M), char(tg(j)), fps, nTracked, nRaw, ratio, V(hit).name}; %#ok<AGROW>
        end
    end

    R = cell2table(rows, 'VariableNames', ...
        {'model','trialTag','fps','trackedFrames','rawFrames','raw_over_tracked','video'});
    disp(R)

    ok = isfinite(R.raw_over_tracked);
    fprintf('\n===== raw / tracked frame ratio by model =====\n');
    if any(ok)
        disp(varfun(@(x)[numel(x) min(x) median(x) max(x)], R(ok,:), ...
             'InputVariables','raw_over_tracked','GroupingVariables','model'))
        fprintf('ratio ~1  -> Stage A kept the whole video (capture was short)\n');
        fprintf('ratio >>1 -> Stage A truncated the track (reprocessing will help)\n');
    else
        fprintf('no videos matched by name; inspect the "video" column\n');
    end

    assignin('base','SRC',R);
    fprintf('\ntable written to base workspace as SRC\n');
end
