function info = parse_trial_name(videoName)
% PARSE_TRIAL_NAME  Infer trial metadata from a video/trial name.
%
%   '25mm_T01'          -> 25 mm drop, trial 1, container 'full'
%   '25mm_T01_shallow'  -> 25 mm drop, trial 1, container 'shallow'
%
%   Returns a struct:
%       ok            logical, true if the name matched
%       dropHeight_mm numeric
%       trialNum      numeric
%       container     'full' | 'shallow'
%       heightLabel   e.g. '25mm'
%       trialParent   e.g. '25mm_T01'        (folder grouping full+shallow)
%       trialTag      e.g. '25mm_T01_full'   (self-describing file prefix)
%       raw           original input
%       reason        why it failed to parse (when ok == false)

    info = struct('ok',false,'dropHeight_mm',NaN,'trialNum',NaN, ...
                  'container','','heightLabel','','trialParent','', ...
                  'trialTag','','raw',char(videoName),'reason','');

    [~, name] = fileparts(char(videoName));   % strip extension if present
    name = strtrim(name);

    % <NN>mm _ T<NN> [ _shallow | _full ]   (tolerant of _ - or space)
    tok = regexpi(name, ...
        '^(\d+)\s*mm[_\-\s]*T(\d+)(?:[_\-\s]*(shallow|full))?$', ...
        'tokens', 'once');

    if isempty(tok)
        info.reason = sprintf(...
            'name "%s" does not match <NN>mm_T<NN>[_shallow|_full]', name);
        return;
    end

    info.dropHeight_mm = str2double(tok{1});
    info.trialNum      = str2double(tok{2});
    if numel(tok) >= 3 && ~isempty(tok{3})
        info.container = lower(tok{3});
    else
        info.container = 'full';
    end

    info.heightLabel = sprintf('%dmm', round(info.dropHeight_mm));
    info.trialParent = sprintf('%dmm_T%02d', round(info.dropHeight_mm), info.trialNum);
    info.trialTag    = sprintf('%s_%s', info.trialParent, info.container);
    info.ok          = true;
end
