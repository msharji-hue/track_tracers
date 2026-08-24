function R = extract_d0_static(varargin)
%EXTRACT_D0_STATIC  Settled penetration of the zero-drop trials, measured.
%
%   READ-ONLY with respect to the results tree: it reads *_tracks.mat and
%   writes only a CSV and a QA figure to OutDir.
%
%   WHY THIS EXISTS. The h = 0 trials are quarantined in load_kinematics_set
%   because kd_kinematics is a DYNAMIC-IMPACT pipeline: its impact/stop
%   detector has no impact to find on a quasi-static release, so the
%   d_final_meas_cm it reports for these trials is an artefact and must not be
%   used. This is the dedicated extraction the loader's BACKLOG note asks for.
%   Nothing here reads, repairs or reuses those quarantined values.
%
%   WHAT IS MEASURED. The settled penetration is a DISPLACEMENT between two
%   quiet states of the same trace, not a depth read off an event:
%       d0_static = median(z over the settled plateau)
%                 - median(z over the pre-release plateau)
%   z(t) is the composition-invariant rod displacement (rod_displacement),
%   positive into the bed. Both plateaus are found by the same quiet-window
%   test, so neither end depends on detecting a release.
%
%   USAGE
%       R = extract_d0_static
%       R = extract_d0_static('Root', 'D:\ME_GRANULAB\JerboaImpact')
%       R = extract_d0_static('PlateauTolCm', 0.03, 'Save', false)
%
%   OPTIONS (name-value)
%       'Root'          results root  (default D:\ME_GRANULAB\JerboaImpact)
%       'OutDir'        output folder (default <Root>/03_RESULTS/_exports)
%       'PlateauTolCm'  rolling-SD ceiling for "quiet", cm      (default 0.02)
%       'MinBaselineMs' shortest pre-release plateau, ms        (default 50)
%       'MinSettleMs'   shortest settled plateau, ms            (default 100)
%       'SdWindowMs'    rolling-SD window, ms                   (default 25)
%       'VReleaseCmS'   a real release must exceed this |dz/dt| (default 1)
%       'Save'          write the CSV and the QA figure         (default true)
%       'Show'          display the QA figure                   (default true)
%
%   OUTPUTS (in OutDir; <stamp> = yyyymmdd_HHMMSS)
%       d0_static_<stamp>.csv      one row per trial
%       d0_static_qa_<stamp>.png   z(t) per trial, both plateaus shaded
%
%   QA FLAGS (semicolon-separated; empty means clean)
%       NOBASELINE  no quiet window of MinBaselineMs before the release
%       NOSETTLE    no quiet window of MinSettleMs at the end
%       NORELEASE   peak |dz/dt| never exceeds VReleaseCmS: nothing moved
%       NEGATIVE    d0 < -0.02 cm, i.e. the foot ended ABOVE where it started
%   Flagged trials are reported and written, and excluded from the summary
%   statistics. They are not repaired here: a flag is a review decision.
%
%   Base MATLAB only -- movstd, median, gradient, exportgraphics.

% ── options ──────────────────────────────────────────────────────────────
opt.Root          = 'D:\ME_GRANULAB\JerboaImpact';
opt.OutDir        = '';
opt.PlateauTolCm  = 0.02;
opt.MinBaselineMs = 50;
opt.MinSettleMs   = 100;
opt.SdWindowMs    = 25;
opt.VReleaseCmS   = 1;
opt.Save          = true;
opt.Show          = true;
optNames = fieldnames(opt);
for i = 1:2:numel(varargin)
    j = find(strcmpi(optNames, varargin{i}), 1);
    if isempty(j)
        error('extract_d0_static:unknownOption', ...
              'Unknown option "%s". Valid: %s', string(varargin{i}), strjoin(optNames', ', '));
    end
    opt.(optNames{j}) = varargin{i+1};
end
if isempty(opt.OutDir)
    opt.OutDir = fullfile(opt.Root, '03_RESULTS', '_exports');
end

D0_DYNAMIC_CM = 0.862;   % Step 3 pure-2/3 multiplicative scale, for comparison ONLY

stamp = datestr(now, 'yyyymmdd_HHMMSS');

fprintf('\n=== extract_d0_static ===\n');
fprintf('  Root        : %s\n', opt.Root);
fprintf('  OutDir      : %s\n', opt.OutDir);
fprintf('  plateau tol : %.3f cm over a %g ms rolling SD\n', opt.PlateauTolCm, opt.SdWindowMs);
fprintf('  min windows : baseline %g ms, settled %g ms\n', opt.MinBaselineMs, opt.MinSettleMs);

% ── 1. locate the zero-drop trials ───────────────────────────────────────
d = dir(fullfile(opt.Root, '03_RESULTS', '**', '*_tracks.mat'));
keep = startsWith({d.name}, '0mm_');
d = d(keep);
if isempty(d)
    error('extract_d0_static:noTrials', ...
          'No 0mm_*_tracks.mat under %s. Zero-drop trials are located by tag.', ...
          fullfile(opt.Root, '03_RESULTS'));
end
n = numel(d);
tag   = strings(n,1);
model = strings(n,1);
for i = 1:n
    tag(i)   = string(regexprep(d(i).name, '_tracks\.mat$', ''));
    model(i) = local_model_from_tag(tag(i));
end
fprintf('\n--- trials ---\n');
fprintf('  zero-drop trials found : %d  (Tight %d, Default %d, Wide %d)\n', ...
        n, sum(model=="Tight"), sum(model=="Default"), sum(model=="Wide"));

% ── 2. per trial: z(t), plateaus, d0 ─────────────────────────────────────
fps      = nan(n,1);   nFrames  = nan(n,1);
d0       = nan(n,1);   baseline = nan(n,1);  settled  = nan(n,1);
baseSd   = nan(n,1);   settSd   = nan(n,1);
vpeak    = nan(n,1);   tSettle  = nan(n,1);
flags    = strings(n,1);
Z        = cell(n,1);  TT = cell(n,1);  BW = cell(n,1);  SW = cell(n,1);

fprintf('\n--- extraction ---\n');
for i = 1:n
    p = fullfile(d(i).folder, d(i).name);
    S = load(p, 'meta', 'tracks', 'calib');

    % fps: the stored rate unless it is outside the plausible band, which is
    % the repair_fps corruption resolve_fps exists to undo.
    f = NaN;
    if isfield(S.tracks,'fps'), f = S.tracks.fps; end
    if ~isfinite(f) || f < 1000 || f > 6000
        f = resolve_fps(p, S.meta, S.tracks);
    end
    if ~isfinite(f)
        flags(i) = "NOFPS";
        fprintf('  %-26s NOFPS: no plausible frame rate, skipped\n', tag(i));
        continue
    end
    fps(i) = f;  dt = 1/f;

    % Composition-invariant rod displacement, positive into the bed. NOT a
    % single marker: markers drop out mid-trial and a positional mean would
    % step when the averaged set changes.
    z = rod_displacement(S.tracks.trackedX, S.tracks.trackedY, ...
                         S.calib.lineA, S.calib.lineB, S.calib.lineC, S.calib.mmPerPx);
    z = z(:).';
    nFrames(i) = numel(z);
    t = (0:numel(z)-1) * dt;
    Z{i} = z;  TT{i} = t;

    % Quiet windows: rolling SD below tolerance, contiguous, long enough.
    wSd   = max(3, round(opt.SdWindowMs*1e-3 / dt));
    sd    = movstd(z, wSd);
    quiet = isfinite(sd) & sd < opt.PlateauTolCm;
    runs  = local_runs(quiet);

    fl = strings(0,1);

    % baseline = FIRST long-enough quiet run; settled = LAST one.
    minB = max(2, round(opt.MinBaselineMs*1e-3 / dt));
    minS = max(2, round(opt.MinSettleMs  *1e-3 / dt));
    bi = [];  si = [];
    if ~isempty(runs)
        len = runs(:,2) - runs(:,1) + 1;
        bi = find(len >= minB, 1, 'first');
        si = find(len >= minS, 1, 'last');
    end

    if isempty(bi)
        fl(end+1,1) = "NOBASELINE"; %#ok<AGROW>
    else
        BW{i} = runs(bi,:);
        baseline(i) = median(z(runs(bi,1):runs(bi,2)));
        baseSd(i)   = std(z(runs(bi,1):runs(bi,2)));
    end
    if isempty(si)
        fl(end+1,1) = "NOSETTLE"; %#ok<AGROW>
    else
        SW{i} = runs(si,:);
        settled(i) = median(z(runs(si,1):runs(si,2)));
        settSd(i)  = std(z(runs(si,1):runs(si,2)));
    end
    % A settled window that IS the baseline window means nothing ever moved
    % away and came back; d0 would be identically zero by construction.
    if ~isempty(bi) && ~isempty(si) && bi == si
        fl(end+1,1) = "NOSETTLE"; %#ok<AGROW>
    end

    v = gradient(z, dt);
    vpeak(i) = max(abs(v));
    if vpeak(i) < opt.VReleaseCmS
        fl(end+1,1) = "NORELEASE"; %#ok<AGROW>
    end

    if ~isempty(bi) && ~isempty(si) && si > bi
        d0(i) = settled(i) - baseline(i);
        % Release onset: first sample after the baseline plateau that exceeds
        % the release speed. Settle time runs from there to the start of the
        % settled plateau.
        k = find(abs(v(runs(bi,2):end)) > opt.VReleaseCmS, 1, 'first');
        if ~isempty(k)
            tSettle(i) = (runs(si,1) - (runs(bi,2)+k-1)) * dt;
        end
    end
    if isfinite(d0(i)) && d0(i) < -0.02
        fl(end+1,1) = "NEGATIVE"; %#ok<AGROW>
    end

    flags(i) = strjoin(fl', ';');
    fprintf('  %-26s %-8s fps %6.0f  d0 = %6s cm  vpeak %6.2f cm/s  %s\n', ...
            tag(i), model(i), fps(i), local_num(d0(i)), vpeak(i), flags(i));
end

% ── 3. summary over unflagged trials ─────────────────────────────────────
ok = flags == "" & isfinite(d0);
fprintf('\n--- d0_static summary (unflagged trials only) ---\n');
fprintf('  %-8s %3s %8s %8s %8s %8s\n', 'model', 'n', 'mean', 'SD', 'min', 'max');
for m = ["Tight" "Default" "Wide"]
    s = ok & model == m;
    if any(s)
        fprintf('  %-8s %3d %8.4f %8.4f %8.4f %8.4f\n', ...
                m, sum(s), mean(d0(s)), std(d0(s)), min(d0(s)), max(d0(s)));
    else
        fprintf('  %-8s %3d %8s %8s %8s %8s\n', m, 0, '-', '-', '-', '-');
    end
end
fprintf('  %-8s %3d %8.4f %8.4f %8.4f %8.4f\n', 'POOLED', sum(ok), ...
        mean(d0(ok)), std(d0(ok)), min(d0(ok)), max(d0(ok)));
fprintf('  pooled d0_static = %.4f +/- %.4f cm (n = %d, %d flagged and excluded)\n', ...
        mean(d0(ok)), std(d0(ok)), sum(ok), n - sum(ok));
fprintf(['  comparison: dynamically inferred scale %.3f cm (Step 3 pure-2/3 fit) ' ...
         'vs measured %.4f cm -- difference %+.4f cm; comparison only, nothing is fitted here\n'], ...
        D0_DYNAMIC_CM, mean(d0(ok)), mean(d0(ok)) - D0_DYNAMIC_CM);

fi = find(flags ~= "");
fprintf('\n--- QA flags ---\n');
if isempty(fi)
    fprintf('  none\n');
else
    for i = fi(:)'
        fprintf('  %-26s %-8s %s\n', tag(i), model(i), flags(i));
    end
end

% ── 4. outputs ───────────────────────────────────────────────────────────
written = strings(0,1);
if opt.Save
    if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
    csvPath = fullfile(opt.OutDir, sprintf('d0_static_%s.csv', stamp));
    fid = fopen(csvPath, 'w');
    if fid < 0
        error('extract_d0_static:csvFailed', 'Could not write %s', csvPath);
    end
    fprintf(fid, '# extract_d0_static.m on %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# d0_static_cm = median(z over settled plateau) - median(z over pre-release plateau)\n');
    fprintf(fid, '# quiet = %g ms rolling SD < %g cm; baseline >= %g ms, settled >= %g ms\n', ...
            opt.SdWindowMs, opt.PlateauTolCm, opt.MinBaselineMs, opt.MinSettleMs);
    fprintf(fid, '# the quarantined dynamic d_final_meas_cm is NOT read or used\n');
    fprintf(fid, ['trialTag,model,fps,nFrames,d0_static_cm,baseline_cm,settled_cm,' ...
                  'baseline_sd,settled_sd,vpeak_cm_s,t_settle_s,qa_flags\n']);
    for i = 1:n
        fprintf(fid, '%s,%s,%.1f,%d,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                tag(i), model(i), fps(i), nFrames(i), local_num(d0(i)), ...
                local_num(baseline(i)), local_num(settled(i)), local_num(baseSd(i)), ...
                local_num(settSd(i)), local_num(vpeak(i)), local_num(tSettle(i)), flags(i));
    end
    fclose(fid);
    written = [written; string(csvPath)];
end

fig = local_figure_qa(tag, Z, TT, BW, SW, d0, flags, opt);
if opt.Save
    pngPath = fullfile(opt.OutDir, sprintf('d0_static_qa_%s.png', stamp));
    exportgraphics(fig, pngPath, 'Resolution', 150);
    written = [written; string(pngPath)];
    fprintf('\n--- written ---\n');
    fprintf('  %s\n', written);
end
if ~opt.Show, close(fig); fig = gobjects(1); end

fprintf('\n');
R = struct();
R.trialTag      = tag;
R.model         = model;
R.fps           = fps;
R.nFrames       = nFrames;
R.d0_static_cm  = d0;
R.baseline_cm   = baseline;
R.settled_cm    = settled;
R.baseline_sd   = baseSd;
R.settled_sd    = settSd;
R.vpeak_cm_s    = vpeak;
R.t_settle_s    = tSettle;
R.qa_flags      = flags;
R.ok            = ok;
R.written       = written;
R.fig           = fig;
end

% ═════════════════════════════════════════════════════════════════════════
function m = local_model_from_tag(tg)
%LOCAL_MODEL_FROM_TAG  Geometry from the tag suffix, as the pipeline names it.
    if endsWith(tg, '_tight')
        m = "Tight";
    elseif endsWith(tg, '_wide')
        m = "Wide";
    else
        m = "Default";
    end
end

function runs = local_runs(mask)
%LOCAL_RUNS  Start/end indices of every contiguous true run in a logical row.
%   Returns an M-by-2 matrix, one run per row. Empty when there are none.
    mask = mask(:).';
    dm = diff([false, mask, false]);
    s = find(dm ==  1);
    e = find(dm == -1) - 1;
    runs = [s(:), e(:)];
end

function s = local_num(v)
%LOCAL_NUM  CSV/console number that stays empty rather than printing NaN.
    if isfinite(v), s = sprintf('%.4f', v); else, s = ''; end
end

function fig = local_figure_qa(tag, Z, TT, BW, SW, d0, flags, opt)
%LOCAL_FIGURE_QA  z(t) per trial with both plateau windows shaded.
%   Small multiples so every trial is inspected, not just the summary: the
%   plateau choice IS the measurement and has to be visible.
    n  = numel(tag);
    nc = ceil(sqrt(n));  nr = ceil(n/nc);
    fig = figure('Color','w','Units','inches','Position',[1 1 2.0*nc 1.5*nr], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl = tiledlayout(fig, nr, nc, 'Padding','compact', 'TileSpacing','compact');
    for i = 1:n
        ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
        if isempty(Z{i})
            text(ax, 0.5, 0.5, 'no trace', 'Units','normalized', ...
                 'HorizontalAlignment','center', 'FontSize', 7);
            set(ax,'XTick',[],'YTick',[]);
            title(ax, char(tag(i)), 'FontSize', 6, 'FontWeight','normal', ...
                  'Interpreter','none');
            continue
        end
        z = Z{i};  t = TT{i};
        yl = [min(z) max(z)];
        if diff(yl) < 1e-6, yl = yl + [-0.05 0.05]; end
        yl = yl + 0.08*diff(yl)*[-1 1];
        % Plateau windows first, so the trace draws over them.
        W = {BW{i}, SW{i}};
        for w = 1:2
            if ~isempty(W{w})
                xw = t(W{w});
                patch(ax, [xw(1) xw(2) xw(2) xw(1)], [yl(1) yl(1) yl(2) yl(2)], ...
                      [0.85 0.90 0.97], 'EdgeColor','none');
            end
        end
        plot(ax, t, z, 'k-', 'LineWidth', 0.5);
        ylim(ax, yl); xlim(ax, [t(1) t(end)]);
        set(ax, 'FontSize', 6, 'LineWidth', 0.4, 'Layer', 'top');
        ttl = sprintf('%s  d_0=%s', tag(i), local_num(d0(i)));
        if flags(i) ~= "", ttl = sprintf('%s  [%s]', ttl, flags(i)); end
        title(ax, ttl, 'FontSize', 6, 'FontWeight','normal', 'Interpreter','none');
    end
    xlabel(tl, 't (s)', 'FontSize', 8);
    ylabel(tl, 'z (cm, positive into the bed)', 'FontSize', 8);
    title(tl, sprintf(['Zero-drop settled penetration: pre-release and settled plateaus ' ...
                       'shaded (tol %.3f cm)'], opt.PlateauTolCm), 'FontSize', 9);
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
