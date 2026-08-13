function [fps, src] = resolve_fps(tracksPath, meta, tracks)
% RESOLVE_FPS  Best per-trial fps, robust to repair_fps corruption.
%
%   repair_fps can overwrite a correct meta.fps_true with a bad median (e.g.
%   2799 -> 600). The tracking scalars CSV (<tag>_scalars.csv) is written BEFORE
%   repair and is never rewritten, so it holds the authoritative track-time rate.
%
%   Cascade, first value inside the plausible band [LO,HI] wins:
%     1) <tag>_scalars.csv  'fps'   (pre-repair truth; fixes the 600 corruption)
%     2) meta.fps_true              (repaired value; used only if CSV implausible,
%                                     i.e. the genuine case repair_fps was for)
%     3) tracks.fps
%   Returns NaN (and src='none') if nothing plausible is found — callers then
%   skip the trial gracefully instead of crashing on an absurd dt.

    LO = 1000; HI = 6000;                 % Hz; rejects the 600 ffmpeg cap.
                                          %  Upper bound raised from 5000 to cover
                                          %  the three-model runs (3300-5200 fps).
                                          %  Real rig rate ~2700-3100, so this
                                          %  band excludes the corruption while
                                          %  leaving ample margin around the truth.
    fps = NaN; src = 'none';

    % 1) pre-repair scalars CSV
    csvp = regexprep(tracksPath, '_tracks\.mat$', '_scalars.csv');
    v = read_csv_fps(csvp);
    if inband(v,LO,HI), fps = v; src = 'scalars_csv'; return; end

    % 2) repaired meta
    if isstruct(meta) && isfield(meta,'fps_true') && inband(meta.fps_true,LO,HI)
        fps = meta.fps_true; src = 'meta.fps_true'; return;
    end

    % 3) tracks struct
    if nargin>=3 && isstruct(tracks) && isfield(tracks,'fps') && inband(tracks.fps,LO,HI)
        fps = tracks.fps; src = 'tracks.fps'; return;
    end
end

function tf = inband(x, lo, hi)
    tf = isscalar(x) && isfinite(x) && x>=lo && x<=hi;
end

function v = read_csv_fps(csvp)
    v = NaN;
    if exist(csvp,'file')~=2, return; end
    try
        fid = fopen(csvp,'r'); h = fgetl(fid); d = fgetl(fid); fclose(fid);
        if ischar(h) && ischar(d)
            H = strsplit(h, ','); D = strsplit(d, ',');
            j = find(strcmpi(strtrim(H),'fps'), 1);
            if ~isempty(j) && j<=numel(D), v = str2double(strtrim(D{j})); end
        end
    catch
        v = NaN;
    end
end
