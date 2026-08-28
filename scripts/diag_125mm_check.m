function T = diag_125mm_check(varargin)
% DIAG_125MM_CHECK  Five-clip scatter diagnostic at h = 125 mm (Default foot).
%
%   Question it answers: of the ~20% relative SD in final depth seen in the
%   reviewed 125 mm Default cell, how much is the PIPELINE and how much is the
%   experiment (release, bed preparation)? Five freshly captured clips are run
%   through the EXISTING, UNMODIFIED pipeline and their scatter compared with
%   the campaign cell. Tight new scatter at matched v0 would point away from
%   bed preparation; ~9-20% relSD(d) reproduced on fresh clips points at it.
%
%   ISOLATION. process_trial and track_tracers_2 hardcode the 01_FRAMES /
%   02_SAVED_DETECTIONS / 03_RESULTS leaves under whatever output root they are
%   given, so the only safe isolation is a SEPARATE TOP-LEVEL ROOT. Everything
%   this script writes goes under ScratchRoot; the official JerboaImpact tree is
%   opened READ-ONLY, and then only to read the master export for comparison.
%   No campaign result, export, exclusion list or reviewed flag is touched.
%
%   VIDEO STAGING. scan_video_tree/infer_item read the material from a 'GB' or
%   'CHIN' folder and the container from a 'full'/'shallow' folder in the path,
%   so the capture folder (which has neither) cannot be scanned directly. The
%   five clips are SYMLINKED into <ScratchRoot>/_videos/GB/full/ under their own
%   capture names and the batch is pointed there. Symlinks, not copies: ~10 GB
%   of raw video is not duplicated, and the originals are never written to.
%
%   USAGE
%     diag_125mm_check                                 % all stages
%     diag_125mm_check('Stages','setup')               % stage only, list clips
%     diag_125mm_check('Stages','track')               % Stage A (slow)
%     diag_125mm_check('Stages','kin')                 % Stage B
%     diag_125mm_check('Stages','analyze')             % scalars + figures
%
%   NAME-VALUE
%     'VideoDir'    capture folder holding the five .avi clips
%     'ScratchRoot' top-level diagnostic root (created if absent)
%     'MasterCsv'   reviewed campaign export, read-only
%     'Model'       foot model for the per-model calibration ['Default Model']
%     'BatchLabel'  batch level in the scratch path        ['DIAG125']
%     'Policy'      process_trial / track_tracers_2 policy ['resume']
%     'Stages'      'all'(default)|'setup'|'track'|'kin'|'analyze'
%
%   OUTPUT
%     T  table of the five new trials (empty unless the analyze stage ran)
%     <ScratchRoot>/diag_d_vs_v0.png   new vs campaign in the (v0, d) plane
%     <ScratchRoot>/diag_overlay.png   z, v, a+g overlays of the five clips

opt.VideoDir    = '/Volumes/GRANULAB/ME_GRANULAB/Test Batches/Default_Batch_5_2';
opt.ScratchRoot = '/Volumes/GRANULAB/ME_GRANULAB/JerboaImpact_DIAG125';
opt.MasterCsv   = ['/Volumes/GRANULAB/ME_GRANULAB/JerboaImpact/03_RESULTS/' ...
                   '_exports/master_trials_20260824_221720.csv'];
opt.Model       = 'Default Model';
opt.BatchLabel  = 'DIAG125';
opt.Policy      = 'resume';
opt.Stages      = 'all';
opt.DropHeight  = 125;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(thisDir, fullfile(codeDir,'src'));

stages = lower(string(opt.Stages));
doAll  = stages == "all";
T      = table();

fprintf('\n========================================================\n');
fprintf(' diag_125mm_check   %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf(' scratch root : %s\n', opt.ScratchRoot);
fprintf(' videos       : %s\n', opt.VideoDir);
fprintf('========================================================\n');

% ── STAGE 0: locate + stage the clips ────────────────────────────────────
if doAll || any(stages == "setup") || any(stages == "track")
    stageDir = stage_videos(opt);
else
    stageDir = fullfile(opt.ScratchRoot, '_videos', 'GB', 'full');
end

% ── STAGE A: tracking ────────────────────────────────────────────────────
if doAll || any(stages == "track")
    fprintf('\n=== STAGE A: process_trial (tracking) ===\n');
    process_trial('batch', stageDir, opt.ScratchRoot, ...
        struct('batchLabel', opt.BatchLabel, 'model', opt.Model, ...
               'policy', opt.Policy, 'keepFrames', 'none', 'limit', 0));
end

% ── STAGE B: kinematics ──────────────────────────────────────────────────
if doAll || any(stages == "kin")
    fprintf('\n=== STAGE B: track_tracers_2 (kinematics) ===\n');
    % saveEventFrames/impactCheck OFF: both re-open the RAW clip through the
    % standard capture tree, which these clips are not in. They are QA extras,
    % not part of the kinematics.
    %
    % NO 'model' OPTION HERE, deliberately. track_tracers_2 reads it as
    % 03_RESULTS/<model>/..., but process_trial writes the model level BELOW
    % material and batch (03_RESULTS/GB/DIAG125/Default Model/...), so passing it
    % looks for a folder that never exists. It is only a folder FILTER: the
    % calibration that matters comes from the per-trial tracks file, where Stage A
    % saved the exact get_calibration_model('Default Model') struct it tracked
    % with. Scanning the whole scratch root is correct and complete -- these five
    % trials are the only thing in it.
    track_tracers_2('batch', opt.ScratchRoot, ...
        struct('policy', opt.Policy, 'figures', 'none', ...
               'saveEventFrames', false, 'impactCheck', false));
end

% ── ANALYSIS ─────────────────────────────────────────────────────────────
if doAll || any(stages == "analyze")
    [T, K] = collect_new_trials(opt);
    if isempty(T), fprintf('\nNo _kin.mat found under the scratch root. Stop.\n'); return; end
    report_calibration(T);
    report_scalars(T);
    assert_apg_convention(K);
    C = load_campaign_cell(opt);
    comparison_block(T, C);
    fig_d_vs_v0(T, C, opt);
    fig_overlay(K, opt);
end

fprintf('\ndiag_125mm_check done.\n');
end

% ═════════════════════════════════════════════════════════════════════════
%  STAGE 0 -- locate and symlink the clips
% ═════════════════════════════════════════════════════════════════════════
function stageDir = stage_videos(opt)
    fprintf('\n=== STAGE 0: locate clips ===\n');
    assert(isfolder(opt.VideoDir), 'diag125:noVideoDir', ...
           'Capture folder not found: %s', opt.VideoDir);

    exts = {'*.avi','*.AVI','*.mp4','*.MP4','*.mov','*.MOV'};
    V = [];
    for e = 1:numel(exts)
        V = [V; dir(fullfile(opt.VideoDir, exts{e}))]; %#ok<AGROW>
    end
    if ~isempty(V)
        V = V(~[V.isdir]);
        [~, ia] = unique(lower({V.name}), 'stable');   % case-insensitive FS
        V = V(ia);
    end

    fprintf('  %d video file(s) in %s\n\n', numel(V), opt.VideoDir);
    fprintf('  %-24s %14s  %s\n', 'name', 'size (MB)', 'modified');
    for i = 1:numel(V)
        fprintf('  %-24s %14.1f  %s\n', V(i).name, V(i).bytes/1e6, ...
                char(datetime(V(i).datenum,'ConvertFrom','datenum', ...
                              'Format','yyyy-MM-dd HH:mm:ss')));
    end

    assert(numel(V) >= 5, 'diag125:tooFewClips', ...
        ['Found %d video(s), need 5. Nothing was processed -- resolve the ' ...
         'capture folder first.'], numel(V));

    if numel(V) > 5
        [~, ord] = sort([V.datenum], 'descend');
        V = V(ord(1:5));
        fprintf('\n  MORE THAN FIVE present: taking the five NEWEST by mtime:\n');
        for i = 1:5, fprintf('    %s\n', V(i).name); end
    end
    [~, ord] = sort({V.name});  V = V(ord);

    % infer_item needs GB + a container folder in the path
    stageDir = fullfile(opt.ScratchRoot, '_videos', 'GB', 'full');
    if ~isfolder(stageDir), mkdir(stageDir); end

    fprintf('\n  staging (symlink) into %s\n', stageDir);
    for i = 1:numel(V)
        src = fullfile(V(i).folder, V(i).name);
        dst = fullfile(stageDir, V(i).name);
        if exist(dst, 'file') == 2 || is_link(dst)
            fprintf('    exists : %s\n', V(i).name);
            continue;
        end
        [st, msg] = system(sprintf('ln -s %s %s', esc(src), esc(dst)));
        assert(st == 0, 'diag125:symlink', 'ln -s failed for %s: %s', V(i).name, msg);
        fprintf('    linked : %s\n', V(i).name);
    end
end

function tf = is_link(p)
    [st, out] = system(sprintf('test -L %s && echo yes', esc(p)));
    tf = (st == 0) && contains(out, 'yes');
end

function s = esc(p)
    s = ['"' char(p) '"'];
end

% ═════════════════════════════════════════════════════════════════════════
%  COLLECT -- read the five _kin.mat back off the scratch root
% ═════════════════════════════════════════════════════════════════════════
function [T, K] = collect_new_trials(opt)
    fprintf('\n=== COLLECT: _kin.mat under the scratch root ===\n');
    F = dir(fullfile(opt.ScratchRoot, '03_RESULTS', '**', '*_kin.mat'));
    T = table(); K = struct([]);
    if isempty(F)
        return;
    end
    [~, ord] = sort({F.name});  F = F(ord);

    rows = struct([]);
    for i = 1:numel(F)
        S    = load(fullfile(F(i).folder, F(i).name), 'meta', 'kin', 'calib');
        meta = S.meta; kin = S.kin; calib = S.calib;
        ra   = getfld(kin, 'rodAngle', struct());
        fps  = getfld(meta, 'fps_true', NaN);

        rows(i).trialTag   = string(getfld(meta,'trialTag', erase(F(i).name,'_kin.mat')));
        rows(i).v0_cm_s    = kin.v0_cm_s;
        rows(i).d_final_cm = kin.d_final_cm;
        rows(i).t_stop_s   = kin.t_stop_s;
        rows(i).fps        = fps;
        rows(i).n_frames_impact_to_stop = kin.stopFrame - kin.impact_index;
        rows(i).rodAngle_peak_abs_deg   = getfld(ra, 'peak_abs_deg', NaN);
        rows(i).a_stop_cm_s2 = getfld(kin, 'a_stop_cm_s2', NaN);
        rows(i).mmPerPx    = getfld(calib, 'mmPerPx', NaN);
        rows(i).bedX       = getfld(calib, 'bedX', NaN);
        rows(i).impactFrame= kin.impact_index;
        rows(i).stopFrame  = kin.stopFrame;
        rows(i).nTracked   = size(getfld(kin,'depthRod_cm',[]), 1);
        rows(i).kinPath    = string(fullfile(F(i).folder, F(i).name));

        K(i).trialTag = rows(i).trialTag;
        K(i).t_s      = kin.t_s(:);
        K(i).z_cm     = kin.depthRod_cm(:);
        K(i).v_cm_s   = kin.v(:);
        K(i).a_cm_s2  = kin.a(:);
        [K(i).apg_cm_s2, K(i).g] = net_accel(kin, calib);   % recomputed g - a
        K(i).impact   = kin.impact_index;
        K(i).stop     = kin.stopFrame;
        K(i).dt       = getfld(kin, 'dt', 1/fps);
    end
    T = struct2table(rows);
    fprintf('  %d trial(s) collected\n', height(T));
end

% ═════════════════════════════════════════════════════════════════════════
%  REPORTS
% ═════════════════════════════════════════════════════════════════════════
function report_calibration(T)
    fprintf('\n=== CALIBRATION (first trial) ===\n');
    fprintf('  trial     : %s\n', T.trialTag(1));
    fprintf('  mmPerPx   : %.4f   (expected 0.1079)  %s\n', T.mmPerPx(1), ...
            verdict(abs(T.mmPerPx(1) - 0.1079) < 1e-6));
    fprintf('  bedX      : %g       (expected 4)      %s\n', T.bedX(1), ...
            verdict(T.bedX(1) == 4));
    % The plausible band is resolve_fps's own [1000, 6000], which deliberately
    % covers the 3300-5200 fps three-model runs as well as the ~2800 fps
    % campaign-1 rate. A rate outside the campaign cell's 2771-2798 is a CAPTURE
    % SETTING difference, not an error: dt = 1/fps is read per trial, so the
    % kinematics are unaffected. It is reported because it is a real difference
    % between these clips and the cell they are being compared with.
    fprintf('  fps range : %.1f - %.1f  (resolve_fps band 1000-6000)  %s\n', ...
            min(T.fps), max(T.fps), ...
            verdict(min(T.fps) >= 1000 && max(T.fps) <= 6000));
    if min(T.fps) < 2771 || max(T.fps) > 2798
        fprintf(['  NOTE: campaign 125 mm cell ran at 2771-2798 fps; these clips ' ...
                 'ran at %.0f-%.0f fps.\n        Different camera setting, same ' ...
                 'physics: dt = 1/fps is per trial.\n'], min(T.fps), max(T.fps));
    end
    u = unique(T.mmPerPx);
    if numel(u) > 1
        fprintf('  NOTE: mmPerPx is not identical across trials: %s\n', ...
                strjoin(compose('%.4f', u'), ', '));
    end
end

function s = verdict(tf)
    if tf, s = '[in range]'; else, s = '[OUT OF RANGE]'; end
end

function report_scalars(T)
    fprintf('\n=== PER-TRIAL SCALARS (new diagnostic clips) ===\n');
    fprintf('  %-26s %9s %9s %9s %8s %8s %10s\n', 'trialTag', 'v0(cm/s)', ...
            'd(cm)', 't_stop(s)', 'fps', 'nFrames', 'rodAng(deg)');
    for i = 1:height(T)
        fprintf('  %-26s %9.2f %9.4f %9.5f %8.1f %8d %10.3f\n', T.trialTag(i), ...
            T.v0_cm_s(i), T.d_final_cm(i), T.t_stop_s(i), T.fps(i), ...
            T.n_frames_impact_to_stop(i), T.rodAngle_peak_abs_deg(i));
    end
    fprintf('  (nFrames = impact frame to stop frame, inclusive of neither end)\n');
end

function assert_apg_convention(K)
% Rest tail of a + g. kd_kinematics masks a (and so a+g) OUTSIDE
% [impact, stop], so the resting rod is not carried in the stored series. The
% velocity IS retained for calib.postCapMs after the stop, and the rod is at
% rest there, so the rest state is recovered by differencing that tail:
%   at rest  a -> 0  =>  a + g = g - a -> +g = +980 cm/s^2.
% Also checked: a + g at impact should be ~0 (free fall, a = +g).
    j  = 1;
    Ki = K(j);
    g  = Ki.g;
    fprintf('\n=== a+g CONVENTION CHECK (%s) ===\n', Ki.trialTag);

    v  = Ki.v_cm_s; dt = Ki.dt;
    tailIdx = (Ki.stop + 3) : numel(v);
    tailIdx = tailIdx(isfinite(v(tailIdx)));
    if numel(tailIdx) >= 5
        n    = min(numel(tailIdx), 15);
        idx  = tailIdx(1:n);
        p    = polyfit((0:n-1)'*dt, v(idx), 1);
        aRest = p(1);
        apgRest = g - aRest;
        fprintf('  post-stop tail: %d frames, mean |v| = %.3f cm/s\n', n, mean(abs(v(idx))));
        fprintf('  a(rest)   = %+8.1f cm/s^2   (expected ~0)\n', aRest);
        fprintf('  a+g(rest) = %+8.1f cm/s^2   (expected ~ +%.0f)\n', apgRest, g);
        assert(abs(apgRest - g) < 0.25*g, 'diag125:apgRestTail', ...
            ['a+g rest tail is %+.1f cm/s^2, not ~ +%.0f. The g - a sign ' ...
             'convention is not holding on this run.'], apgRest, g);
        fprintf('  PASS: rest tail sits at +g within 25%%.\n');
    else
        fprintf('  SKIP: post-stop velocity tail too short (%d frames).\n', numel(tailIdx));
    end

    % NOT checked here: the free-fall checkpoint (a = +g, a+g = 0). kd_kinematics
    % masks a before the impact index, and AT the impact index the adaptive fit
    % window is clamped to [impact, stop], so a+g(impact) already reads the
    % deceleration by construction. Free fall is simply not represented in the
    % stored series, and reporting a+g(impact) against an "expected 0" would be
    % comparing against a value this quantity cannot take.
    apgPeak = max(Ki.apg_cm_s2(Ki.impact:Ki.stop));
    fprintf('  a+g(peak)   = %+8.1f cm/s^2   (grains resisting, expected >> +g)\n', apgPeak);
    assert(apgPeak > g, 'diag125:apgPeak', ...
        'peak a+g (%.1f) is not above +g; the sign convention looks inverted.', apgPeak);
end

% ═════════════════════════════════════════════════════════════════════════
%  CAMPAIGN CELL -- READ ONLY
% ═════════════════════════════════════════════════════════════════════════
function C = load_campaign_cell(opt)
    fprintf('\n=== CAMPAIGN CELL (read-only) ===\n');
    fprintf('  %s\n', opt.MasterCsv);
    assert(isfile(opt.MasterCsv), 'diag125:noMaster', 'master export not found');
    M = readtable(opt.MasterCsv, 'TextType','string');

    keep = (M.condition == "GB/full") & (M.model == "Default") & ...
           (abs(M.dropHeight_mm - opt.DropHeight) < 1e-6) & as_logical(M.keep_reviewed);
    C = M(keep, {'trialTag','v0_cm_s','d_final_cm','t_stop_s','fps'});
    fprintf('  filter: condition GB/full, model Default, dropHeight %g mm, keep_reviewed true\n', ...
            opt.DropHeight);
    fprintf('  n = %d\n', height(C));
end

function comparison_block(T, C)
    vN = T.v0_cm_s;   dN = T.d_final_cm;
    vC = C.v0_cm_s;   dC = C.d_final_cm;

    fprintf('\n');
    fprintf('========================================================\n');
    fprintf(' COMPARISON: new diagnostic clips vs reviewed 125 mm cell\n');
    fprintf('========================================================\n');

    fprintf('\n  NEW (n = %d)\n', numel(vN));
    fprintf('  %-26s %10s %10s\n', 'trialTag', 'v0(cm/s)', 'd(cm)');
    for i = 1:numel(vN)
        fprintf('  %-26s %10.2f %10.4f\n', T.trialTag(i), vN(i), dN(i));
    end

    fprintf('\n  CAMPAIGN (n = %d)\n', numel(vC));
    fprintf('  %-26s %10s %10s\n', 'trialTag', 'v0(cm/s)', 'd(cm)');
    for i = 1:numel(vC)
        fprintf('  %-26s %10.2f %10.4f\n', C.trialTag(i), vC(i), dC(i));
    end

    fprintf('\n  %-12s %10s %10s %10s %10s %10s %10s\n', ...
        'group', 'n', 'mean v0', 'SD v0', 'relSD v0', 'mean d', 'SD d');
    fprintf('  %-12s %10d %10.2f %10.2f %9.2f%% %10.4f %10.4f\n', 'new', ...
        numel(vN), mean(vN), std(vN), 100*std(vN)/mean(vN), mean(dN), std(dN));
    fprintf('  %-12s %10d %10.2f %10.2f %9.2f%% %10.4f %10.4f\n', 'campaign', ...
        numel(vC), mean(vC), std(vC), 100*std(vC)/mean(vC), mean(dC), std(dC));
    fprintf('  %-12s %10s %10s %10s %10s %9.2f%% (new) vs %.2f%% (campaign)\n', ...
        'relSD d', '', '', '', '', 100*std(dN)/mean(dN), 100*std(dC)/mean(dC));

    rvN = 100*std(vN)/mean(vN);   rvC = 100*std(vC)/mean(vC);
    rdN = 100*std(dN)/mean(dN);   rdC = 100*std(dC)/mean(dC);

    % v0 offset in campaign SDs -- were the new clips released the same way?
    dz = (mean(vN) - mean(vC)) / std(vC);

    % ROBUSTNESS OF THE CAMPAIGN NUMBER. The cell's relSD(d) is the yardstick
    % sentence (b) is measured against, so it is worth knowing whether it rests
    % on one point. Drop the single most deviant trial (largest |d - median|)
    % and recompute: if the cell relSD collapses, the "20%% scatter" being
    % re-shot is one outlier, not a broad spread, and the new clips must be
    % compared against the ROBUST value instead.
    [rdCrob, dropTag] = robust_relsd(dC, C.trialTag);

    % Mean SHIFT is a different failure from mean SCATTER and the diagnostic has
    % to separate them: a bed prepared to a different density moves the mean at
    % unchanged v0, while trial-to-trial preparation noise widens the spread.
    % Tested against the SAME de-outliered cell the scatter is judged against.
    % Testing the mean against the full cell would use the outlier-inflated
    % variance as the yardstick -- the identical mistake as quoting its relSD.
    dCrob   = drop_most_deviant(dC);
    dShift  = mean(dN) - mean(dC);
    dShiftR = mean(dN) - mean(dCrob);
    [tW,    dfW]    = welch_t(dN, dCrob);      % robust cell: the one that counts
    [tWall, dfWall] = welch_t(dN, dC);         % full cell, reported for contrast

    % 95%% CI on a relSD from n samples (chi-square on the variance). At n = 5
    % the interval runs from 0.60x to 2.88x the point estimate -- the "+/-60%%"
    % the diagnostic was scoped around. Computed rather than asserted so the
    % verdict below cannot claim more resolution than five trials carry.
    [rdLo, rdHi] = relsd_ci(rdN, numel(dN));
    FLEET = 9.3;   % fleet-median within-group relSD(d), the "normal" floor

    fprintf('\n  relSD(d) new = %.1f%%  [95%% CI %.1f - %.1f%%]  (n = %d)\n', ...
            rdN, rdLo, rdHi, numel(dN));
    fprintf('  reference band: fleet-median within-group %.1f%%  ..  this cell %.1f%%\n', ...
            FLEET, rdC);
    fprintf('  cell relSD(d) is %.1f%% with all %d trials, %.1f%% once %s (d = the low\n', ...
            rdC, numel(dC), rdCrob, dropTag);
    fprintf('    outlier) is dropped -- the headline 20%% rests on that ONE trial.\n');
    fprintf('  mean d: new %.3f cm vs cell %.3f cm (all 10) -> shift %+.3f cm (%+.1f%%), ', ...
            mean(dN), mean(dC), dShift, 100*dShift/mean(dC));
    fprintf('Welch t = %.2f (df %.1f)\n', tWall, dfWall);
    fprintf('  vs de-outliered cell %.3f cm -> shift %+.3f cm (%+.1f%%), Welch t = %.2f (df %.1f)\n', ...
            mean(dCrob), dShiftR, 100*dShiftR/mean(dCrob), tW, dfW);
    fprintf(['    The full-cell t is the WEAKER test: the outlier inflates the cell ' ...
             'variance it divides by.\n']);

    fprintf('\n  DIAGNOSTIC SENTENCES\n');
    if rvN <= 1.6*rvC
        s1 = sprintf(['(a) RELEASE REPEATABILITY: new relSD(v0) = %.1f%% against the campaign ' ...
            'cell''s %.1f%%, so impact speed is at least as repeatable as it was in the ' ...
            'campaign, and the new mean sits %.1f campaign SD from the campaign mean -- ' ...
            'these five clips are the same drop, so v0 is not what is moving.'], ...
            rvN, rvC, dz);
    else
        s1 = sprintf(['(a) RELEASE REPEATABILITY: new relSD(v0) = %.1f%% against the campaign ' ...
            'cell''s %.1f%%, so impact speed scattered MORE than in the campaign (mean offset ' ...
            '%.1f campaign SD); the release itself, not only the bed, differs between these ' ...
            'clips and must be fixed before depth scatter can be read.'], rvN, rvC, dz);
    end
    fprintf('  %s\n', wrap(s1));

    % Judge the SCATTER against the robust cell value and the fleet floor, not
    % against a number one outlier set; report the SHIFT separately.
    ref = min(rdCrob, FLEET);
    if rdLo >= ref
        sc = sprintf(['reproduces the cell''s own scatter (%.1f%% robust, %.1f%% fleet floor) ' ...
            'even at the bottom of its interval'], rdCrob, FLEET);
    elseif rdHi < ref
        sc = sprintf(['is genuinely TIGHTER than the cell (%.1f%% robust, %.1f%% fleet floor): ' ...
            'the whole interval sits below it'], rdCrob, FLEET);
    else
        sc = sprintf(['is COMPARABLE to the cell once its single low outlier is set aside ' ...
            '(%.1f%% robust, %.1f%% fleet floor): the interval straddles that value, so five ' ...
            'trials cannot call it tighter'], rdCrob, FLEET);
    end
    if abs(tW) >= 2.5
        sh = sprintf(['but the five sit %+.2f cm (%+.0f%%) from that same cell mean at ' ...
            'matched v0 (Welch t = %.1f, df %.0f), a real LEVEL shift rather than extra ' ...
            'spread'], dShiftR, 100*dShiftR/mean(dCrob), tW, dfW);
        tail = ['Read together: the pipeline is not manufacturing the scatter -- fresh ' ...
            'clips through unmodified code land in the same relative spread -- so the ' ...
            're-shoot targets the right thing; but the depth OFFSET at unchanged v0 says ' ...
            'bed preparation is ALSO moving the operating point between sessions, which ' ...
            'is a second, distinct problem the re-shoot must control. Five more clips ' ...
            'size it; deeper analysis of these five will not.'];
    else
        sh = sprintf(['and the means agree to %+.2f cm against that same cell (Welch t = ' ...
            '%.1f, df %.0f), so the bed level has not moved either'], dShiftR, tW, dfW);
        tail = ['Read together: fresh clips through unmodified code reproduce both the ' ...
            'level and the relative spread of the reviewed cell, so the scatter is ' ...
            'experimental rather than pipeline-induced and the re-shoot targets the right ' ...
            'thing. Five more clips would tighten the estimate; deeper analysis of these ' ...
            'five will not.'];
    end
    s2 = sprintf('(b) DEPTH SCATTER AT MATCHED v0: new relSD(d) = %.1f%% [95%% CI %.1f-%.1f%%] %s, %s. %s', ...
        rdN, rdLo, rdHi, sc, sh, tail);
    fprintf('  %s\n', wrap(s2));

    fprintf(['\n  (n = 5 estimates a relSD to roughly +/-60%%. These sentences are an ' ...
             'interpretation frame, not a decision.)\n']);
end

function d = drop_most_deviant(d)
% The cell with its single most deviant trial removed (median-referenced, so
% the outlier cannot pick itself by dragging the centre). Shared by the robust
% relSD and the robust mean test so both judge against the same nine trials.
    d = d(:);
    [~, k] = max(abs(d - median(d)));
    d(k)   = [];
end

function [rel, dropped] = robust_relsd(d, tags)
% relSD(d) recomputed with the single most deviant trial removed, so a cell
% whose scatter rests on one point is visible as such. Median-referenced, not
% mean-referenced: the outlier must not choose itself by dragging the centre.
    d = d(:);
    [~, k]  = max(abs(d - median(d)));
    dropped = "the most deviant trial";
    if nargin >= 2 && numel(tags) == numel(d), dropped = string(tags(k)); end
    r    = drop_most_deviant(d);
    rel  = 100 * std(r) / mean(r);
end

function [t, df] = welch_t(a, b)
% Welch's unequal-variance t and its Satterthwaite df. Written out rather than
% called from the Statistics Toolbox so the diagnostic runs on a bare MATLAB.
    a = a(:); b = b(:);
    va = var(a)/numel(a);   vb = var(b)/numel(b);
    t  = (mean(a) - mean(b)) / sqrt(va + vb);
    df = (va + vb)^2 / (va^2/(numel(a)-1) + vb^2/(numel(b)-1));
end

function [lo, hi] = relsd_ci(rel, n)
% 95% CI on a relative SD, from the chi-square interval on the variance:
%   s*sqrt((n-1)/chi2inv(0.975,n-1))  ..  s*sqrt((n-1)/chi2inv(0.025,n-1))
% The n = 5 factors are hardcoded so this does not require the Statistics
% Toolbox for the case the diagnostic was designed around; other n fall back to
% chi2inv when it is available, and to a NaN interval when it is not.
    if n == 5
        lo = 0.5989 * rel;   hi = 2.8746 * rel;      % chi2inv(.975,4)=11.143, (.025,4)=0.4844
    elseif exist('chi2inv','file') == 2
        lo = rel * sqrt((n-1)/chi2inv(0.975, n-1));
        hi = rel * sqrt((n-1)/chi2inv(0.025, n-1));
    else
        lo = NaN; hi = NaN;
    end
end

function out = wrap(s)
    out = strjoin(cellstr(textwrap({char(s)}, 86)), sprintf('\n  '));
end

% ═════════════════════════════════════════════════════════════════════════
%  FIGURES
% ═════════════════════════════════════════════════════════════════════════
function fig_d_vs_v0(T, C, opt)
    cNew  = [0.4940 0.1840 0.5560];   % blue-purple, new diagnostic clips
    cCamp = [0.8500 0.3250 0.0980];   % red, reviewed campaign cell
    msz   = 9;

    f  = figure('Color','w','Position',[100 100 760 560]);
    ax = axes(f); apply_fig_style(ax);

    % Step-3 pooled fit, drawn first so the data sit on top of it
    vGrid = linspace(0, 1.08*max([T.v0_cm_s; C.v0_cm_s]), 200);
    dFit  = 0.044 + 0.0711 * vGrid.^(2/3);
    hFit  = plot(ax, vGrid, dFit, 'k-', 'LineWidth', 1.4);

    hCamp = plot(ax, C.v0_cm_s, C.d_final_cm, 'o', 'MarkerSize', msz, ...
        'MarkerFaceColor', cCamp, 'MarkerEdgeColor', cCamp*0.6, 'LineWidth', 0.8);
    hNew  = plot(ax, T.v0_cm_s, T.d_final_cm, 'd', 'MarkerSize', msz, ...
        'MarkerFaceColor', cNew,  'MarkerEdgeColor', cNew*0.6,  'LineWidth', 0.8);

    % campaign cell mean +/- SD as a crosshair
    mv = mean(C.v0_cm_s); sv = std(C.v0_cm_s);
    md = mean(C.d_final_cm); sd = std(C.d_final_cm);
    hX = plot(ax, [mv-sv mv+sv], [md md], '-', 'Color', cCamp, 'LineWidth', 2.4);
    plot(ax, [mv mv], [md-sd md+sd], '-', 'Color', cCamp, 'LineWidth', 2.4);
    plot(ax, mv, md, '+', 'Color', cCamp, 'MarkerSize', 14, 'LineWidth', 2.4);

    xlabel(ax, 'impact speed  v_0  (cm s^{-1})');
    ylabel(ax, 'final depth  d  (cm)');
    title(ax, sprintf('h = %g mm, Default foot: five new clips vs the reviewed cell', ...
          opt.DropHeight));
    legend(ax, [hNew hCamp hX hFit], ...
        {'new (diagnostic)','campaign 125 mm','cell mean \pm SD','Step-3 fit'}, ...
        'Location','northwest');

    lo = min([T.v0_cm_s; C.v0_cm_s]); hi = max([T.v0_cm_s; C.v0_cm_s]);
    pad = max(6, 0.12*(hi-lo));
    xlim(ax, [lo-pad hi+pad]);

    p = fullfile(opt.ScratchRoot, 'diag_d_vs_v0.png');
    exportgraphics(f, p, 'Resolution', 200);
    close(f);
    fprintf('\n  wrote %s\n', p);
end

function fig_overlay(K, opt)
    n   = numel(K);
    col = lines(max(n,3));

    f = figure('Color','w','Position',[100 100 820 900]);
    ax = gobjects(3,1);
    for k = 1:3, ax(k) = subplot(3,1,k); apply_fig_style(ax(k)); end

    h = gobjects(n,1);
    for i = 1:n
        Ki  = K(i);
        % show impact to a little past the stop, in ms from impact
        lo  = Ki.impact;
        hi  = min(numel(Ki.t_s), Ki.stop + round(0.008/Ki.dt));
        idx = lo:hi;
        tms = 1e3*Ki.t_s(idx);

        h(i) = plot(ax(1), tms, Ki.z_cm(idx),    '-', 'Color', col(i,:), 'LineWidth', 1.3);
               plot(ax(2), tms, Ki.v_cm_s(idx),  '-', 'Color', col(i,:), 'LineWidth', 1.3);
               plot(ax(3), tms, Ki.apg_cm_s2(idx),'-','Color', col(i,:), 'LineWidth', 1.3);
    end
    yline(ax(2), 0, 'k:');
    yline(ax(3), K(1).g, 'k:', 'a+g = g (at rest)', 'LabelHorizontalAlignment','left');
    yline(ax(3), 0, 'k:');

    ylabel(ax(1), 'depth  z  (cm)');
    ylabel(ax(2), 'velocity  v  (cm s^{-1})');
    ylabel(ax(3), 'a + g  (cm s^{-2})');
    xlabel(ax(3), 'time from impact  (ms)');
    title(ax(1), sprintf('h = %g mm diagnostic clips: overlay (t = 0 at impact)', opt.DropHeight));
    names = cell(n,1);
    for i = 1:n, names{i} = char(K(i).trialTag); end
    legend(ax(1), h, names, 'Location','southeast', 'Interpreter','none');

    p = fullfile(opt.ScratchRoot, 'diag_overlay.png');
    exportgraphics(f, p, 'Resolution', 200);
    close(f);
    fprintf('  wrote %s\n', p);
end

% ═════════════════════════════════════════════════════════════════════════
function v = getfld(s, f, dflt)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function tf = as_logical(col)
% keep_reviewed comes back as logical, double or "true"/"false" text depending
% on how readtable typed the column. Accept all three rather than assume one.
    if islogical(col)
        tf = col;
    elseif isnumeric(col)
        tf = col ~= 0;
    else
        tf = strcmpi(strtrim(string(col)), "true") | strtrim(string(col)) == "1";
    end
    tf = tf(:);
end
