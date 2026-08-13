function meta = load_meta_checked(tracksPath, quiet)
% LOAD_META_CHECKED  Load meta from a *_tracks.mat and attach the corrected
% drop height, warning loudly when a known label problem applies.
%
%   *** USE THIS INSTEAD OF load(tp,'meta') IN ANALYSIS CODE. ***
%
%   Filenames, folders, trialTags and meta.dropHeight_mm are deliberately NOT
%   renamed on disk. That means a plain load() silently returns the WRONG height
%   for GB/shallow. This wrapper is the safety net: it adds
%
%       meta.dropHeight_true_mm   physical height (use this for all physics)
%       meta.heightCorrected      true if the label was reversed
%       meta.heightNote           human-readable explanation
%
%   and warns once per MATLAB session so the issue cannot pass unnoticed.
%
%   meta = load_meta_checked(tracksPath)        % warns once per session
%   meta = load_meta_checked(tracksPath, true)  % suppress the session warning

    persistent warned          % must be declared before any executable line
    if nargin < 2, quiet = false; end
    S = load(tracksPath, 'meta');
    meta = S.meta;

    cond = string(getf(meta,'material','')) + "/" + string(getf(meta,'container',''));
    hlab = getf(meta,'dropHeight_mm', NaN);

    [htrue, corrected, note] = true_drop_height(hlab, cond);
    meta.dropHeight_true_mm = htrue;
    meta.heightCorrected    = corrected;
    meta.heightNote         = note;

    if corrected && ~quiet
        if isempty(warned)
            warned = true;
            warning('load_meta_checked:reversedLabels', ...
              ['\n' repmat('=',1,70) '\n' ...
               ' GB/shallow DROP-HEIGHT LABELS ARE REVERSED ON DISK.\n' ...
               ' Filenames/trialTags keep the ORIGINAL (wrong) label by design.\n' ...
               ' Use meta.dropHeight_true_mm for anything physical.\n' ...
               ' e.g. %s : labelled %g mm  ->  ACTUAL %g mm\n' ...
               ' Map lives in src/true_drop_height.m (single source of truth).\n' ...
               repmat('=',1,70)], getf(meta,'trialTag','?'), hlab, htrue);
        end
    end
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
