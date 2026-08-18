function videoPath = find_raw_video(rawRoot, meta)
%FIND_RAW_VIDEO  Locate the raw .avi for a trial. Single source of truth.
%
%   videoPath = find_raw_video(rawRoot, meta)
%     rawRoot  ancestor of the capture tree; searched recursively
%     meta     needs heightLabel, trialNum, material, container
%
%   MATCHING IS EXACT on the filename stem AND requires the material and
%   container folders in the path. A substring test on the height label is not
%   enough: '25mm' is contained in '125mm' and '325mm', which is how an earlier
%   version opened CHIN/as_poured/125mm_T04.avi when asked for 25mm_T04_dense.
%
%   The ORIGINAL capture is preferred over any '_transcoded' copy, because the
%   transcoded AVIs carry the 600 fps muxer artefact that resolve_fps exists to
%   reject.
%
%   Errors if nothing matches; the message names the stem and the folders that
%   were required, so a wrong rawRoot is obvious.

    V = dir(fullfile(rawRoot, '**', '*.avi'));
    V = V(~[V.isdir]);
    if isempty(V)
        error('find_raw_video:noVideos', ...
              'No .avi found anywhere under %s. Check rawRoot (JERBOA_RAW_ROOT).', rawRoot);
    end
    % Windows lists *.avi and *.AVI twice; unique collapses that.
    paths = unique(string(fullfile({V.folder}', {V.name}')));

    stem = sprintf('%s_T%02d', meta.heightLabel, meta.trialNum);
    [~, stems] = arrayfun(@(p) fileparts(p), paths, 'UniformOutput', false);
    stems = string(stems);

    matDir  = string(filesep) + string(meta.material)  + string(filesep);
    contDir = string(filesep) + string(meta.container) + string(filesep);

    hit = strcmpi(stems, stem) & contains(paths, matDir,  'IgnoreCase', true) ...
                               & contains(paths, contDir, 'IgnoreCase', true);
    if ~any(hit)
        error('find_raw_video:noMatch', ...
            'No raw video with stem "%s" under a %s%s path in %s', ...
            stem, matDir, contDir, rawRoot);
    end

    cand = paths(hit);
    orig = cand(~contains(cand, "transcoded", 'IgnoreCase', true));
    if ~isempty(orig), cand = orig; end
    videoPath = char(cand(1));
end
