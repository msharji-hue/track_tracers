function it = infer_item(fp)
% INFER_ITEM  Infer all trial metadata for ONE video, from its path + name.
%
%   material  : nearest 'GB' or 'CHIN' folder in the path
%   container : nearest 'full' or 'shallow' folder in the path
%               (folder wins; falls back to the filename suffix; else 'full')
%   drop/trial: parsed from the filename (parse_trial_name)
%
%   Returns a struct: fullpath,name,folder,material,container,dropHeight_mm,
%   trialNum,heightLabel,trialParent,trialTag,ok,reason.

    fp = char(fp);
    [folder, nm, ext] = fileparts(fp);
    name = [nm ext];

    parts     = local_split(folder);
    material  = local_pick(parts, {'GB','CHIN'});
    container = local_pick(parts, {'full','shallow','as_poured','dense'});
    info      = parse_trial_name(name);

    it = struct('fullpath',fp, 'name',name, 'folder',folder, ...
        'material',upper(material), 'container','', ...
        'dropHeight_mm',NaN, 'trialNum',NaN, 'heightLabel','', ...
        'trialParent','', 'trialTag','', 'ok',false, 'reason','');

    % container precedence: folder -> filename -> 'full'
    if ~isempty(container)
        it.container = lower(container);
    elseif info.ok
        it.container = info.container;
    else
        it.container = 'full';
    end

    if isempty(material)
        it.reason = 'material (GB/CHIN) not found in path';
        return;
    end
    if ~info.ok
        it.reason = info.reason;
        return;
    end

    it.dropHeight_mm = info.dropHeight_mm;
    it.trialNum      = info.trialNum;
    it.heightLabel   = info.heightLabel;
    it.trialParent   = info.trialParent;
    it.trialTag      = sprintf('%s_%s', info.trialParent, it.container);
    it.ok            = true;
end

% ── helpers ───────────────────────────────────────────────────────────────
function parts = local_split(p)
    parts = regexp(char(p), '[\\/]+', 'split');
    parts = parts(~cellfun('isempty', parts));
end

function tok = local_pick(parts, set)
    tok = '';
    for i = numel(parts):-1:1            % search deepest folder first
        for j = 1:numel(set)
            if strcmpi(parts{i}, set{j})
                tok = set{j};
                return;
            end
        end
    end
end
