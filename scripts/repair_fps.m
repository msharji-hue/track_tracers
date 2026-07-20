function repair_fps(outRoot, rawRoot, varargin)
% REPAIR_FPS  Restore the true frame rate on trials processed from transcoded
% video. METADATA ONLY — no re-tracking, no re-export, no re-detection.
%
% ROOT CAUSE
%   ffmpeg's AVI muxer will not write a ~2778 fps rate: it silently falls back
%   to 600 fps. (Verified: -r, -video_track_timescale and -time_base all still
%   produce 600/1 in the AVI header.) Every FRAME is preserved, so detection and
%   tracking are unaffected — but VideoReader.FrameRate returned 600, and
%   process_one_trial stored that into meta.fps_true / tracks.fps and the
%   _scalars.csv fps column for every trial processed from the transcoded tree.
%   Re-transcoding CANNOT fix this; only the stored metadata can be corrected.
%
% WHY IT MATTERS
%   Downstream kinematics use dt = 1/fps. At 600 instead of ~2778, every
%   velocity is wrong by ~4.6x and every acceleration by ~21x — silently,
%   because the tracks themselves look perfect.
%
% HOW THE TRUE fps IS RECOVERED, in order of preference
%   1. ffprobe on the ORIGINAL raw video (the camera wrote a correct rate;
%      ffmpeg can read these files even though VideoReader could not).
%   2. Median fps of sibling trials in the same (material, condition, height)
%      cell that already carry a good rate.
%   3. Median of the same (material, condition).
%   4. Global median.
%   Anything not resolved is left untouched and reported as UNRESOLVED.
%
%   The provenance of every repaired value is written to
%   <outRoot>/03_RESULTS/_batch_logs/fps_repair_<timestamp>.csv so the choice
%   is auditable and can be stated in the methods section.
%
% USAGE
%   repair_fps(outRoot, rawRoot)                    % DRY RUN (default) — reports only
%   repair_fps(outRoot, rawRoot, 'dryRun', false)   % apply the repair
%
%   outRoot : .../ME_GRANULAB/JerboaImpact
%   rawRoot : .../ME_GRANULAB/Test Batches/Batch 5   (the ORIGINAL tree, not _transcoded)
%
% NOTE  The raw videos may be Dropbox online-only placeholders. ffprobe will
%       trigger a download per file. If that is slow, make the raw tree
%       available offline first, or rely on the sibling-median fallback.

    p = inputParser;
    addParameter(p,'dryRun',       true);
    addParameter(p,'suspectBelow', 1000);       % true rates are ~2700-3100
    addParameter(p,'ffprobe',      '/opt/homebrew/bin/ffprobe');
    parse(p, varargin{:});
    o = p.Results;

    fprintf('\n=== REPAIR FPS (%s) ===\n', ternary(o.dryRun,'DRY RUN — nothing written','APPLYING'));
    fprintf('  outRoot : %s\n  rawRoot : %s\n\n', outRoot, rawRoot);

    % -- 1) inventory every tracked trial (partial load: meta only, fast) ----
    files = dir(fullfile(outRoot, '03_RESULTS', '**', 'tracks', '*_tracks.mat'));
    if isempty(files)
        error('repair_fps: no _tracks.mat found under %s', fullfile(outRoot,'03_RESULTS'));
    end
    T = struct('path',{},'tag',{},'material',{},'container',{}, ...
               'height',{},'trial',{},'fps',{});
    for i = 1:numel(files)
        pth = fullfile(files(i).folder, files(i).name);
        try
            S = load(pth, 'meta');
        catch
            fprintf('  WARNING: could not read meta from %s — skipped\n', files(i).name);
            continue;
        end
        m = S.meta;
        T(end+1) = struct('path',pth, 'tag',m.trialTag, ...
            'material',m.material, 'container',m.container, ...
            'height',round(m.dropHeight_mm), 'trial',m.trialNum, ...
            'fps',fld(m,'fps_true'));  %#ok<AGROW>
    end

    allFps  = [T.fps];
    isGood  = isfinite(allFps) & allFps >= o.suspectBelow;
    suspect = find(~isGood);

    fprintf('  trials found      : %d\n', numel(T));
    fprintf('  fps OK            : %d  (median %.1f)\n', sum(isGood), median(allFps(isGood)));
    fprintf('  fps SUSPECT (<%d) : %d\n\n', o.suspectBelow, numel(suspect));
    if isempty(suspect)
        fprintf('Nothing to repair.\n'); return;
    end

    % -- 2) resolve a true fps for each suspect ------------------------------
    logRows = cell(numel(suspect),1);
    nFixed = 0; nUnres = 0;
    for k = 1:numel(suspect)
        idx = suspect(k);
        t   = T(idx);

        rawPath = fullfile(rawRoot, t.material, t.container, ...
                           sprintf('%dmm_T%02d.avi', t.height, t.trial));

        [fpsNew, src] = resolve_fps(rawPath, T, isGood, t, o);

        if ~isfinite(fpsNew)
            src = 'UNRESOLVED'; nUnres = nUnres + 1;
            fprintf('  [%3d/%3d] %-28s %6.1f -> %-10s %s\n', ...
                k, numel(suspect), t.tag, t.fps, 'NaN', src);
        else
            fprintf('  [%3d/%3d] %-28s %6.1f -> %8.1f   (%s)\n', ...
                k, numel(suspect), t.tag, t.fps, fpsNew, src);
            if ~o.dryRun
                apply_fps(t.path, fpsNew);
                nFixed = nFixed + 1;
            end
        end

        logRows{k} = sprintf('%s,%s,%s,%d,%d,%.4f,%s,%s,%s', ...
            t.tag, t.material, t.container, t.height, t.trial, ...
            t.fps, num2str(fpsNew,'%.4f'), src, t.path);
    end

    % -- 3) provenance log ---------------------------------------------------
    logDir = fullfile(outRoot, '03_RESULTS', '_batch_logs');
    if ~exist(logDir,'dir'), mkdir(logDir); end
    logPath = fullfile(logDir, sprintf('fps_repair_%s.csv', datestr(now,'yyyymmdd_HHMMSS')));
    fid = fopen(logPath,'w');
    fprintf(fid, 'trialTag,material,condition,dropHeight_mm,trialNum,fps_old,fps_new,fps_source,tracksPath\n');
    for k = 1:numel(logRows), fprintf(fid, '%s\n', logRows{k}); end
    fclose(fid);

    fprintf('\n  repaired   : %d\n', nFixed);
    fprintf('  unresolved : %d\n', nUnres);
    fprintf('  log        : %s\n', logPath);
    if o.dryRun
        fprintf('\n  DRY RUN — nothing was written. Re-run with ''dryRun'',false to apply.\n');
    end
end

% ------------------------------------------------------------------------
function [fps, src] = resolve_fps(rawPath, T, isGood, t, o)
% Preference: ffprobe raw > sibling cell median > condition median > global.
    fps = NaN; src = '';

    % (1) ffprobe the ORIGINAL raw video
    if isfile(rawPath)
        f = probe_fps(o.ffprobe, rawPath);
        if isfinite(f) && f >= o.suspectBelow
            fps = f; src = 'ffprobe_raw'; return;
        end
    end

    % (2) same material + condition + drop height
    sel = isGood & strcmp({T.material}, t.material) & ...
                   strcmp({T.container}, t.container) & ([T.height] == t.height);
    if any(sel), fps = median([T(sel).fps]); src = 'median_cell'; return; end

    % (3) same material + condition
    sel = isGood & strcmp({T.material}, t.material) & strcmp({T.container}, t.container);
    if any(sel), fps = median([T(sel).fps]); src = 'median_condition'; return; end

    % (4) global
    if any(isGood), fps = median([T(isGood).fps]); src = 'median_global'; end
end

function fps = probe_fps(ffprobe, videoPath)
% Read the true frame rate from a video header via ffprobe. Returns NaN on any
% failure. Tries r_frame_rate first, then avg_frame_rate.
    fps = NaN;
    if ~isfile(ffprobe), return; end
    cmd = sprintf(['%s -v error -select_streams v:0 ' ...
                   '-show_entries stream=r_frame_rate,avg_frame_rate ' ...
                   '-of default=noprint_wrappers=1 "%s"'], ffprobe, videoPath);
    [status, out] = system(cmd);
    if status ~= 0, return; end
    keys = {'r_frame_rate','avg_frame_rate'};
    for i = 1:numel(keys)
        tok = regexp(out, [keys{i} '=(\d+)/(\d+)'], 'tokens', 'once');
        if ~isempty(tok)
            den = str2double(tok{2});
            if den > 0
                v = str2double(tok{1}) / den;
                if isfinite(v) && v > 0, fps = v; return; end
            end
        end
    end
end

function apply_fps(matPath, fpsNew)
% Rewrite the trial's saved outputs with the corrected rate, reusing the
% existing save_tracks() so the file layout and CSV schema are unchanged.
% Only the fps fields change; trackedX/trackedY are untouched.
    S = load(matPath);                            % meta, tracks, calib
    S.meta.fps_true = fpsNew;
    S.tracks.fps    = fpsNew;
    resultsDir = fileparts(fileparts(matPath));   % .../<trial>/<condition>
    save_tracks(S.meta, S.tracks, S.calib, resultsDir);
end

function v = fld(s, f)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
