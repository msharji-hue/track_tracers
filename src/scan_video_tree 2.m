function items = scan_video_tree(rootDir)
% SCAN_VIDEO_TREE  Recursively find every video under rootDir and infer
% per-video metadata via infer_item (material/container from folders,
% drop/trial from filename).
%
%   Handles the structure:
%       <root>/GB/full/25mm_T01.mp4
%       <root>/GB/shallow/25mm_T01.mp4
%       <root>/CHIN/full/...
%   and also works if you point it directly at a leaf folder.
%
%   Requires MATLAB R2016b+ for the '**' recursive glob.

    rootDir = char(rootDir);
    exts    = {'*.avi','*.AVI','*.mp4','*.MP4','*.mov','*.MOV'};

    files = [];
    for e = 1:numel(exts)
        files = [files; dir(fullfile(rootDir, '**', exts{e}))]; %#ok<AGROW>
    end

    items = repmat(infer_item('seed'), 0, 1);   % typed empty
    if isempty(files), return; end

    % de-duplicate (case-insensitive filesystems list *.mp4 and *.MP4 twice)
    full = arrayfun(@(f) fullfile(f.folder, f.name), files, 'UniformOutput', false);
    [~, ia] = unique(lower(full), 'stable');
    files   = files(ia);

    for k = 1:numel(files)
        items(k,1) = infer_item(fullfile(files(k).folder, files(k).name)); %#ok<AGROW>
    end

    % stable ordering: material, container, drop height, trial
    keys = arrayfun(@(it) sprintf('%-5s|%-8s|%07.1f|%05.0f', ...
        it.material, it.container, nz(it.dropHeight_mm), nz(it.trialNum)), ...
        items, 'UniformOutput', false);
    [~, order] = sort(keys);
    items = items(order);
end

function x = nz(v)
    if isnan(v), x = 0; else, x = v; end
end
