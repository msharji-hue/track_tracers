function [h_true_mm, wasCorrected, note] = true_drop_height(h_label_mm, condition)
% TRUE_DROP_HEIGHT  Canonical labelled -> actual drop-height correction.
%
%   *** SINGLE SOURCE OF TRUTH. Do not hardcode this map anywhere else. ***
%
%   GB/shallow drop-height labels were RECORDED IN REVERSE ORDER. The labels
%   (and therefore every filename, folder name, trialTag and meta.dropHeight_mm)
%   are left UNCHANGED on disk; this function supplies the physical height.
%
%   Verified by v0: with the labelled heights, measured v0/sqrt(2gh) reached
%   3.9-4.2 for GB/shallow (physically impossible, >1 means faster than free
%   fall) and the measured v0 was nearly CONSTANT (~275, 265, 240 cm/s) across a
%   5x range of labelled heights. After correction every ratio falls to ~1.0-1.1.
%
%   The map is a SWAP and therefore self-inverse (25<->365, 65<->325, 125<->285,
%   165 fixed). Applying it twice returns the original, so it must be applied
%   exactly once, at read time, from this function only.
%
%   [h_true_mm, wasCorrected, note] = true_drop_height(h_label_mm, condition)
%     h_label_mm : labelled height(s), mm (scalar or array)
%     condition  : "GB/shallow", "GB/full", "CHIN/dense", "CHIN/as_poured"
%                  (also accepts a bare container: "shallow", "full", ...)
%
%   Only GB/shallow is altered. All other conditions pass through unchanged.

    LABEL  = [ 25  65 125 165 285 325 365];
    ACTUAL = [365 325 285 165 125  65  25];

    h_label_mm = double(h_label_mm);
    h_true_mm  = h_label_mm;
    n          = numel(h_label_mm);
    wasCorrected = false(size(h_label_mm));

    cs = lower(string(condition));
    if isscalar(cs) && n > 1, cs = repmat(cs, size(h_label_mm)); end

    for i = 1:n
        if ~contains(cs(min(i,numel(cs))), "shallow"), continue; end
        j = find(LABEL == h_label_mm(i), 1);
        if isempty(j)
            warning('true_drop_height:unmappedHeight', ...
                ['GB/shallow labelled height %g mm is not in the correction map ' ...
                 '(%s). Left unchanged - verify manually.'], ...
                h_label_mm(i), mat2str(LABEL));
            continue;
        end
        h_true_mm(i)    = ACTUAL(j);
        wasCorrected(i) = true;
    end

    if nargout >= 3
        if any(wasCorrected)
            note = "GB/shallow: drop-height labels are REVERSED on disk; " + ...
                   "corrected via true_drop_height.m";
        else
            note = "";
        end
    end
end
