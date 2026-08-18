function videoPath = find_raw_video(rawRoot, meta, model)
%FIND_RAW_VIDEO  Locate the raw .avi for a trial. Single source of truth.
%
%   videoPath = find_raw_video(rawRoot, meta)
%   videoPath = find_raw_video(rawRoot, meta, model)
%     rawRoot  ancestor of the capture tree; searched recursively
%     meta     needs heightLabel, trialNum, material, container
%     model    foot model for this trial. Defaults to meta.model. Accepts any
%              spelling: 'Tight Model', 'Tight', 'tight', 'Tight_Model'.
%
%   MATCHING, in order:
%
%   1. EXACT filename stem, AND the material and container folders in the path.
%      A substring test on the height label is not enough: '25mm' is contained
%      in '125mm' and '325mm', which is how an earlier version opened
%      CHIN/as_poured/125mm_T04.avi when asked for 25mm_T04_dense.
%
%   2. THE MODEL FOLDER. The stem alone is NOT unique: the same 25mm_T04.avi
%      exists under Default, Tight and Wide, so stem+material+container matched
%      three different clips and the first was returned silently. A caller then
%      seeks a frame that may not exist in that clip -- crashing on CurrentTime
%      if it is shorter, or, worse, displaying the wrong trial's frame if it is
%      not. So when the capture tree is organised by model, the trial's own
%      model folder is required.
%
%      Whether the tree IS model-organised is decided from the candidates
%      themselves, not assumed: if none of them sits under a recognised model
%      folder, the tree predates the model layout (campaign 1) and the filter is
%      skipped rather than rejecting every candidate.
%
%   3. The ORIGINAL capture is preferred over any '_transcoded' copy, because
%      the transcoded AVIs carry the 600 fps muxer artefact that resolve_fps
%      exists to reject. That is a duplicate of the SAME clip, so resolving it
%      is not ambiguity.
%
%   4. Anything still ambiguous is an ERROR listing every candidate. Picking the
%      first silently is what caused the bug above; a wrong clip that merely
%      looks plausible is worse than a stop.

    if nargin < 3 || isempty(model)
        model = '';
        if isstruct(meta) && isfield(meta,'model'), model = meta.model; end
    end

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

    % ── 2) model folder ──────────────────────────────────────────────────
    KNOWN = ["default","tight","wide"];
    candModel = arrayfun(@(p) local_model_of_path(p, KNOWN), cand);
    treeHasModels = any(candModel ~= "");

    if treeHasModels
        key = local_model_key(model);
        if key == ""
            error('find_raw_video:modelUnknown', ...
                ['The capture tree under %s is organised by model, but trial ' ...
                 '%s carries no model, so the right clip cannot be chosen. ' ...
                 'Candidates:\n%s'], ...
                rawRoot, stem, local_list(cand));
        end
        keep = candModel == key;
        if ~any(keep)
            error('find_raw_video:noModelMatch', ...
                ['No raw video for stem "%s" under a "%s" model folder in %s. ' ...
                 'Found the same stem under: %s\nThe clip for this model is ' ...
                 'missing, or the model folder is named differently.'], ...
                stem, key, rawRoot, strjoin(cellstr(unique(candModel(candModel~=""))), ', '));
        end
        cand = cand(keep);
    end

    % ── 3) original over transcoded (same clip, not ambiguity) ───────────
    orig = cand(~contains(cand, "transcoded", 'IgnoreCase', true));
    if ~isempty(orig), cand = orig; end

    % ── 4) still ambiguous -> stop, never guess ──────────────────────────
    if numel(cand) > 1
        error('find_raw_video:ambiguous', ...
            ['%d raw videos match trial "%s" (model "%s") and nothing ' ...
             'distinguishes them. Refusing to guess -- an earlier version ' ...
             'returned the first, which silently opened another model''s ' ...
             'clip.\n%s\nNarrow rawRoot, or remove the duplicates.'], ...
            numel(cand), stem, local_model_key(model), local_list(cand));
    end

    videoPath = char(cand(1));
end

% ─────────────────────────────────────────────────────────────────────────
function key = local_model_key(model)
%LOCAL_MODEL_KEY  Normalise any model spelling to 'default'|'tight'|'wide'|''.
%   Mirrors get_calibration_model: strip the word "model", spaces and
%   underscores, lowercase. So 'Tight Model', 'tight', 'Tight_Model' and
%   'TIGHT' all resolve to "tight".
    key = lower(strtrim(string(model)));
    key = erase(key, ["model", "_", " "]);
    if key == "", key = ""; end
end

function k = local_model_of_path(p, known)
%LOCAL_MODEL_OF_PATH  Which model folder, if any, this path sits under.
%   Tests whole path COMPONENTS, so a stray substring in a filename cannot be
%   mistaken for a model folder.
    k = "";
    parts = split(string(p), filesep);
    for i = 1:numel(parts)
        c = local_model_key(parts(i));
        if any(c == known), k = c; return; end
    end
end

function s = local_list(cand)
    s = strjoin("    " + string(cand(:)).', newline);
end
