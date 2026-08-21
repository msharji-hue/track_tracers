function tags = get_manual_exclusions()
%GET_MANUAL_EXCLUSIONS  Trials excluded by hand. Single source of truth.
%
%   tags = get_manual_exclusions()
%
%   Returns the trialTags dropped by deliberate review rather than by an
%   automatic rule. load_kinematics_set applies this list by default, so every
%   analysis that goes through the loader excludes the same trials without
%   each script keeping its own copy -- which is how the lists drifted apart
%   before load_default_gb was retired.
%
%   AUTOMATIC vs MANUAL. The automatic rules (NaN depth, GLITCH, depth cut)
%   are conditions the code can test. This list is for judgements it cannot:
%   a velocity profile a person looked at and rejected. Anything that could be
%   expressed as a rule belongs in the rules, not here.
%
%   EVERY ENTRY CARRIES ITS OWN REASON AND DATE. They are identical today
%   because they came from one review, but they will not stay that way: a tag
%   added next month for a different reason must say so on its own line, and a
%   tag removed later should be removable without disturbing any other.
%
%   TO ADD a trial: append a line with the reason and the date of the decision.
%   TO REMOVE one: delete its line. Do not comment it out -- git history is the
%   record of what was once excluded, and a commented-out tag reads as
%   ambiguous forever.
%
%   Tags are unique across models (the _default / _tight / _wide suffix is part
%   of the tag), so one flat list covers every geometry.

tags = [ ...
    "25mm_T04_full_default"    % stop-fit failure / bad velocity profile, 2026-08-20 review
    "85mm_T09_full_default"    % stop-fit failure / bad velocity profile, 2026-08-20 review
    "105mm_T04_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "285mm_T07_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "285mm_T08_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "305mm_T07_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "305mm_T08_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "325mm_T03_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "325mm_T04_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "345mm_T02_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "345mm_T03_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "345mm_T04_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "345mm_T06_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "345mm_T07_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "345mm_T08_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "345mm_T10_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "345mm_T12_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "345mm_T19_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "365mm_T07_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "365mm_T08_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "365mm_T10_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "365mm_T11_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    "365mm_T12_full_default"   % stop-fit failure / bad velocity profile, 2026-08-20 review
    ];

tags = tags(:);

% A duplicated tag would be counted once by ismember but reported twice in the
% loader's "N of M listed" line, which is exactly the kind of small
% inconsistency that makes a provenance record untrustworthy.
if numel(unique(tags)) ~= numel(tags)
    d = unique(tags(arrayfun(@(t) sum(tags == t) > 1, tags)));
    error('get_manual_exclusions:duplicateTag', ...
        'Tag(s) listed more than once: %s', strjoin(cellstr(d), ', '));
end
end
