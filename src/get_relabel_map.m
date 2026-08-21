function R = get_relabel_map()
%GET_RELABEL_MAP  Read-time corrections to recorded drop heights.
%
%   R = get_relabel_map()
%
%   Returns a two-column table, trialTag -> corrected dropHeight_mm. Every
%   entry is a trial whose recorded height label is believed wrong and whose
%   correct height was established by review.
%
%   READ-TIME, NEVER WRITE-TIME. This is the same house pattern as
%   get_manual_exclusions and the phi refresh in load_kinematics_set: the
%   correction is applied when the analysis loads the data, and NOTHING on
%   disk changes. Raw videos, 01_FRAMES, the tracks, the _kin.mat files and
%   the _kin_scalars.csv files all keep their original labels, including the
%   height baked into the trialTag itself. Consequences that matter:
%
%     * the trialTag is an identifier, not a measurement -- a relabeled trial
%       still reads "25mm_..." forever, and that is deliberate. Renaming files
%       would break provenance back to the capture session and to every
%       diagnostic ever run on that tag.
%     * the correction is reversible by deleting a line here.
%     * anything reading the CSVs directly rather than through
%       load_kinematics_set gets the ORIGINAL label. That is the price of not
%       mutating stored data; load_kinematics_set is the single place the
%       correction is applied, and it records what it did (see the
%       'relabeled' and 'dropHeight_asRecorded' columns it adds).
%
%   EVERY ENTRY CARRIES ITS OWN REASON AND DATE, as in get_manual_exclusions.
%   A relabel is a stronger claim than an exclusion -- it moves a trial INTO
%   another bin, where it then influences that bin's median -- so the evidence
%   for it belongs on the line.
%
%   TO ADD a trial: append a row with the reason and the date of the decision.
%   TO REMOVE one: delete its row. Do not comment it out; git history is the
%   record of what was once relabeled.

% Char vectors, not double-quoted strings: string() on a cell array of
% character vectors is the documented conversion, and this cell exists only to
% keep each tag on the same line as its corrected height and its reason.
%
% trialTag                      corrected height (mm)
E = { ...
    '25mm_T02_full_default',    65   % matches the Default 65 mm bin on v0, d_final and t_stop: d2 = 1.2, margin to next-best bin >= 19. Matcher CSV relabel_candidates.csv, 2026-08-21; trace overlays confirmed.
    '25mm_T03_full_default',    65   % matches the Default 65 mm bin on v0, d_final and t_stop: d2 = 5.7, margin to next-best bin >= 19. Matcher CSV relabel_candidates.csv, 2026-08-21; trace overlays confirmed.
    };

if isempty(E)
    % Deleting the last row must leave a usable empty map, not a table whose
    % two columns disagree about how many rows they have.
    R = table(strings(0,1), zeros(0,1), ...
        'VariableNames', {'trialTag','dropHeight_mm'});
    return
end

R = table(string(E(:,1)), cell2mat(E(:,2)), ...
    'VariableNames', {'trialTag','dropHeight_mm'});

% A tag listed twice would relabel to whichever row ismember happened to find
% first -- a silent, order-dependent choice of drop height. Stop instead.
if numel(unique(R.trialTag)) ~= height(R)
    d = unique(R.trialTag(arrayfun(@(t) sum(R.trialTag == t) > 1, R.trialTag)));
    error('get_relabel_map:duplicateTag', ...
        'Tag(s) listed more than once: %s', strjoin(cellstr(d), ', '));
end

% A non-finite or negative corrected height would propagate into every height
% axis and every fit as a quietly bad bin.
bad = ~isfinite(R.dropHeight_mm) | R.dropHeight_mm < 0;
if any(bad)
    error('get_relabel_map:badHeight', ...
        'Corrected height must be finite and >= 0. Offending tag(s): %s', ...
        strjoin(cellstr(R.trialTag(bad)), ', '));
end
end
