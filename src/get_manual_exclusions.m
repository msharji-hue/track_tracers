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
%   The v0 consistency guard in load_kinematics_set is deliberately NOT a rule
%   for this reason inverted: it can be tested, but what to DO about a
%   violation is a judgement. It warns; the trials review decided to drop
%   because of it are listed here, by name, with a reason.
%
%   EVERY ENTRY CARRIES ITS OWN REASON AND DATE. The list is kept as one
%   block per review decision, concatenated at the end: the block header
%   carries the full reasoning, and each line carries a short reason and the
%   date so it still stands alone if the block is ever split up.
%
%   TO ADD a trial: append a line to the appropriate block, or start a new
%   block if the reason is new. TO REMOVE one: delete its line. Do not comment
%   it out -- git history is the record of what was once excluded, and a
%   commented-out tag reads as ambiguous forever.
%
%   RELABELING IS A DIFFERENT ANSWER. A trial whose height label is wrong but
%   whose data are sound belongs in get_relabel_map, not here: excluding it
%   throws away a good measurement, relabeling keeps it in the bin it actually
%   came from. Excluding is for trials with no bin to go to.
%
%   Tags are unique across models (the _default / _tight / _wide suffix is part
%   of the tag), so one flat list covers every geometry.

% ── 2026-08-20 review ────────────────────────────────────────────────────
% Stop-fit failure or a velocity profile rejected on inspection.
stopFitFailure = [ ...
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

% ── 2026-08-21 integrity review, block A ─────────────────────────────────
% v0/t_stop inconsistent with any ladder height (likely off-ladder drop or
% failed trial).
%
% These five are not individually odd. As a GROUP they share a coherent
% low-v0 / long-t_stop signature that matches no rung of the height ladder --
% and in particular does not match 65 mm, which is where the two relabeled
% 25 mm trials went. A consistent signature that fits nothing is a different
% failure from a noisy trace, and it is why these are excluded rather than
% relabeled: there is no bin to move them to.
offLadder = [ ...
    "25mm_T08_full_default"    % off-ladder / failed trial, 2026-08-21 review
    "25mm_T09_full_default"    % off-ladder / failed trial, 2026-08-21 review
    "25mm_T10_full_default"    % off-ladder / failed trial, 2026-08-21 review
    "25mm_T11_full_default"    % off-ladder / failed trial, 2026-08-21 review
    "25mm_T12_full_default"    % off-ladder / failed trial, 2026-08-21 review
    ];

% ── 2026-08-21 integrity review, block B ─────────────────────────────────
% Measurement failure: v0 inconsistent with the trial's own d_final and/or
% above any free-fall bound.
%
% Two independent checks, either sufficient. The free-fall check is the one
% load_kinematics_set now applies as a warning guard; these are the trials
% review resolved as genuine failures rather than as the documented Tight and
% Wide release-height offsets (see the data notes in README.md). The Tight
% 45 mm trials T03/T05/T06/T09/T10 were considered here and deliberately KEPT:
% they are one population with their clean 45 mm siblings under that offset.
measurementFailure = [ ...
    "145mm_T04_full_default"   % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "165mm_T01_full_default"   % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "185mm_T09_full_default"   % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "285mm_T03_full_default"   % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "325mm_T01_full_default"   % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "45mm_T06_full_default"    % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "45mm_T07_full_default"    % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "85mm_T10_full_default"    % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "145mm_T06_full_tight"     % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "25mm_T02_full_tight"      % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "105mm_T08_full_wide"      % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "225mm_T10_full_wide"      % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "365mm_T04_full_wide"      % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    "85mm_T07_full_wide"       % v0 inconsistent with own d_final / free fall, 2026-08-21 review
    ];

% ── 2026-08-21 integrity review, block C ─────────────────────────────────
% v0 matches the 125 mm class but d_final = 1.71 cm sits outside that bin's
% core, so the two observables disagree about which height this was. No clean
% match anywhere on the ladder; excluded by maintainer decision after overlay
% review rather than relabeled into a bin it only half fits.
noCleanMatch = [ ...
    "105mm_T03_full_default"   % v0/d_final disagree, no clean height match, 2026-08-21 review
    ];

tags = [stopFitFailure; offLadder; measurementFailure; noCleanMatch];

tags = tags(:);

% A duplicated tag would be counted once by ismember but reported twice in the
% loader's "N of M listed" line, which is exactly the kind of small
% inconsistency that makes a provenance record untrustworthy. With the list
% split into blocks this also catches the same trial being excluded twice for
% two different stated reasons -- which would mean one of the reasons is wrong.
if numel(unique(tags)) ~= numel(tags)
    d = unique(tags(arrayfun(@(t) sum(tags == t) > 1, tags)));
    error('get_manual_exclusions:duplicateTag', ...
        'Tag(s) listed more than once: %s', strjoin(cellstr(d), ', '));
end
end
