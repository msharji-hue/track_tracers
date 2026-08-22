function R = relabel_candidates(root, varargin)
%RELABEL_CANDIDATES  Match v0-inconsistent trials against clean height bins.
%   For each flagged trial, score every candidate height bin in the SAME
%   model on z(v0), z(d_final), z(t_stop) vs the bin's clean-trial robust
%   statistics; combined distance d2 = sum(z.^2). READ-ONLY: proposes, never
%   applies. Verdicts:
%     RELABEL    best d2 <= 8, all |z| <= 2.5, margin to 2nd best >= 4
%     AMBIGUOUS  fits best bin but 2nd best within margin
%     NO-MATCH   fails everywhere -> remains excluded
%   Physical gate per candidate: 0.5 <= v0/sqrt(2 g h_cand) <= 1.15.
%   For RELABEL/AMBIGUOUS, overlays the trial's z(t) on the target bin's
%   clean traces (trace-shape evidence).
%
%   R = relabel_candidates(root);            % table + evidence figures
%   Options: 'MakePlots' (true), 'D2Accept' (8), 'ZMax' (2.5), 'Margin' (4)

opt.MakePlots = true; opt.D2Accept = 8; opt.ZMax = 2.5; opt.Margin = 4;
for i=1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
G = 980;  % cm/s^2

flags = [ ...
 "25mm_T02_full_default" "25mm_T03_full_default" "25mm_T08_full_default" ...
 "25mm_T09_full_default" "25mm_T10_full_default" "25mm_T11_full_default" ...
 "25mm_T12_full_default" ...
 "45mm_T03_full_tight" "45mm_T05_full_tight" "45mm_T06_full_tight" ...
 "45mm_T09_full_tight" "45mm_T10_full_tight" ...
 "105mm_T03_full_default" "145mm_T04_full_default" "165mm_T01_full_default" ...
 "185mm_T09_full_default" "285mm_T03_full_default" "325mm_T01_full_default" ...
 "45mm_T06_full_default" "45mm_T07_full_default" "85mm_T10_full_default" ...
 "145mm_T06_full_tight" "25mm_T02_full_tight" ...
 "105mm_T08_full_wide" "225mm_T10_full_wide" "365mm_T04_full_wide" ...
 "85mm_T07_full_wide"]';

% load WITHOUT dropping the flagged trials (they must appear in the table)
K = load_kinematics_set(root, 'exclude', cellstr(get_manual_exclusions()));
K = K(~K.isZeroDrop, :);
isFlag = ismember(K.trialTag, flags);
clean  = K(~isFlag & K.keep, :);

% measurement-scale floors so tight bins don't over-reject
FLOOR = struct('v0', 3.0, 'd', 0.08, 't', 0.0020);   % cm/s, cm, s

rows = {};
for f = flags'
    r = K(K.trialTag == f, :);
    if isempty(r), fprintf('%-26s NOT IN SET\n', f); continue; end
    mdl = r.model(1);
    cand = unique(clean.dropHeight_mm(clean.model == mdl));
    best = []; scores = [];
    for h = cand'
        B = clean(clean.model == mdl & clean.dropHeight_mm == h, :);
        if height(B) < 4, continue; end
        ratio = r.v0_cm_s / sqrt(2*G*h/10);
        if ratio > 1.15 || ratio < 0.5, continue; end
        s = @(x,fl) max(1.4826*mad(x,1), fl);
        z = [ (r.v0_cm_s    - median(B.v0_cm_s))    / s(B.v0_cm_s,    FLOOR.v0); ...
              (r.d_final_cm - median(B.d_final_cm)) / s(B.d_final_cm, FLOOR.d); ...
              (r.t_stop_s   - median(B.t_stop_s))   / s(B.t_stop_s,   FLOOR.t) ];
        scores = [scores; h, z', sum(z.^2)]; %#ok<AGROW>
    end
    if isempty(scores)
        rows{end+1} = {char(f), r.dropHeight_mm, r.v0_cm_s, r.d_final_cm, ...
                       NaN, NaN, NaN, NaN, NaN, NaN, "NO-MATCH (no viable bin)"}; %#ok<AGROW>
        continue
    end
    scores = sortrows(scores, 5);
    b = scores(1,:);  d2_2 = NaN; h2 = NaN;
    if size(scores,1) >= 2, d2_2 = scores(2,5); h2 = scores(2,1); end
    ok  = b(5) <= opt.D2Accept && all(abs(b(2:4)) <= opt.ZMax);
    sep = isnan(d2_2) || (d2_2 - b(5) >= opt.Margin);
    if ok && sep,      v = "RELABEL -> " + b(1) + " mm";
    elseif ok,         v = "AMBIGUOUS (" + b(1) + " vs " + h2 + " mm)";
    else,              v = "NO-MATCH";
    end
    rows{end+1} = {char(f), r.dropHeight_mm, r.v0_cm_s, r.d_final_cm, ...
                   b(1), b(2), b(3), b(4), b(5), d2_2, v}; %#ok<AGROW>
end
R = cell2table(vertcat(rows{:}), 'VariableNames', ...
   {'trialTag','label_mm','v0','d_final','best_mm','z_v0','z_d','z_t', ...
    'd2_best','d2_next','verdict'});
disp(R)
fprintf('\nRELABEL: %d | AMBIGUOUS: %d | NO-MATCH: %d\n', ...
    sum(startsWith(R.verdict,"RELABEL")), sum(startsWith(R.verdict,"AMBIGUOUS")), ...
    sum(startsWith(R.verdict,"NO-MATCH")));
writetable(R, fullfile(root,'03_RESULTS','_batch_logs','relabel_candidates.csv'));

% trace-shape evidence for RELABEL / AMBIGUOUS
if opt.MakePlots
    P = R(~startsWith(R.verdict,"NO-MATCH"), :);
    for i = 1:height(P)
        r  = K(K.trialTag == P.trialTag{i}, :);
        B  = clean(clean.model == r.model(1) & clean.dropHeight_mm == P.best_mm(i), :);
        figure('Name', P.trialTag{i}); hold on; grid on
        for j = 1:height(B)
            kb = load(B.kinPath(j)).kin;                       % adjust field if needed
            plot((kb.t_s - kb.t_s(kb.impact_index))*1000, kb.z, 'Color', [.7 .7 .7]);
        end
        kt = load(r.kinPath(1)).kin;
        plot((kt.t_s - kt.t_s(kt.impact_index))*1000, kt.z, 'r', 'LineWidth', 2);
        title(sprintf('%s over clean %d mm bin (%s)', P.trialTag{i}, ...
              P.best_mm(i), r.model(1)), 'Interpreter','none');
        xlabel('t - t_{impact} (ms)'); ylabel('depth (cm)'); xlim([-5 50]);
    end
end
end
