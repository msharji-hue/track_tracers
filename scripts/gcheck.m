function R = gcheck(root, varargin)
% GCHECK  Calibration/friction GATE driver.
%   Assembles trials from <root>/03_RESULTS/**/*_tracks.mat and runs
%   validate_calibration_gcheck. This is the convenience entry point that the
%   runbook's "Step 3" refers to — validate_calibration_gcheck itself takes an
%   already-assembled trials struct array, not a root path.
%
%   R = gcheck(root)                 % one trial per drop height (quick sweep)
%   R = gcheck(root,'select',{'285mm_T08_full','365mm_T03_dense'})
%   R = gcheck(root,'limit',12)      % first N discovered (sorted)
%   R = gcheck(root,'material','GB') % restrict material (GB|CHIN)
%   R = gcheck(root,'container','full')
%   R = gcheck(root,'massG',65,'gBand',[0.85 1.00])   % passed through
%
%   Fields handed to validate_calibration_gcheck per trial:
%     .trackedX .trackedY .fps .h_cm(=dropHeight_mm/10) .label(=trialTag)

    p = inputParser;
    addParameter(p,'select',{},@iscell);
    addParameter(p,'limit',0,@isnumeric);
    addParameter(p,'material','',@(x)ischar(x)||isstring(x));
    addParameter(p,'container','',@(x)ischar(x)||isstring(x));
    addParameter(p,'massG',65,@(x)isempty(x)||isnumeric(x));
    addParameter(p,'gBand',[0.85 1.00],@isnumeric);
    parse(p,varargin{:});
    o = p.Results;

    thisDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(fileparts(thisDir),'src'));
    calib = get_calibration();

    D = dir(fullfile(root,'03_RESULTS','**','*_tracks.mat'));
    if isempty(D), error('No *_tracks.mat found under %s/03_RESULTS', root); end

    % ── light metadata pass (load meta only) ──────────────────────────────
    M = struct('tag',{},'material',{},'container',{},'h',{},'trial',{},'path',{});
    for k = 1:numel(D)
        tp = fullfile(D(k).folder, D(k).name);
        try s = load(tp,'meta'); m = s.meta; catch, continue; end
        if ~isempty(o.material)  && ~strcmpi(getfld(m,'material',''),  o.material),  continue; end
        if ~isempty(o.container) && ~strcmpi(getfld(m,'container',''), o.container), continue; end
        M(end+1) = struct('tag',getfld(m,'trialTag',''), ...
            'material',getfld(m,'material',''), 'container',getfld(m,'container',''), ...
            'h',getfld(m,'dropHeight_mm',NaN), 'trial',getfld(m,'trialNum',NaN), ...
            'path',tp); %#ok<AGROW>
    end
    if isempty(M), error('No trials matched the material/container filter.'); end

    % stable order: material, container, height, trial
    keys = arrayfun(@(x) sprintf('%s_%s_%08.1f_%04.0f', x.material, x.container, x.h, x.trial), ...
        M, 'UniformOutput', false);
    [~,ord] = sort(keys); M = M(ord);

    % ── selection ─────────────────────────────────────────────────────────
    if ~isempty(o.select)
        M = M(ismember({M.tag}, o.select));
    elseif o.limit > 0
        M = M(1:min(o.limit, numel(M)));
    else
        [~, iu] = unique([M.h], 'stable');   % one trial per distinct drop height
        M = M(iu);
    end
    if isempty(M), error('Selection left no trials.'); end

    % ── load tracks for the chosen few, build the trials struct array ─────
    trials = struct('trackedX',{},'trackedY',{},'fps',{},'h_cm',{},'label',{});
    for j = 1:numel(M)
        s = load(M(j).path,'meta','tracks'); m = s.meta; tr = s.tracks;
        fps = getfld(m,'fps_true',NaN);
        if ~isfinite(fps) && isfield(tr,'fps'), fps = tr.fps; end
        trials(end+1) = struct('trackedX',tr.trackedX,'trackedY',tr.trackedY, ...
            'fps',fps,'h_cm',M(j).h/10,'label',M(j).tag); %#ok<AGROW>
    end
    fprintf('gcheck: %d trials selected from %s\n', numel(trials), root);

    R = validate_calibration_gcheck(trials, calib, 'massG',o.massG, 'gBand',o.gBand);
end

function v = getfld(s, f, dflt)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
