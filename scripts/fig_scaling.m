function R = fig_scaling(varargin)
%FIG_SCALING  The three final manuscript scaling figures.
%
%   READ-ONLY with respect to the results tree: it reads the master trial
%   export and writes only the three figures.
%
%   PRESENTATION WORK ONLY. Every number drawn here is an already-approved
%   result. The coefficients are recomputed with the same least squares the
%   analysis scripts used, purely so a figure cannot drift from the result it
%   claims to show, and each one is ASSERTED against its approved value before
%   anything is drawn. Nothing new is fitted and no selection is revisited.
%
%   STYLE. These deliberately keep the look of step3_fits_working.png: plain
%   MATLAB defaults -- default colour order, default fonts, default line
%   weights, grid on, box on. No house style is applied, because the working
%   figure IS the intended look and any restyling here would silently make the
%   manuscript figures disagree with the figures the analysis was read from.
%   The ONE deliberate visual change is marker SHAPE per geometry, so the
%   three models stay distinguishable in print and to colourblind readers.
%
% ═══════════════════════════════════════════════════════════════════════════
%  CAPTION NOTES  (the display-vs-fit rule; carry this into the captions)
% ═══════════════════════════════════════════════════════════════════════════
%
%   FITS USE N = 524. Every fit, every residual and the RSS(n) profile are
%   computed on the full approved selection: keep_reviewed & ~isZeroDrop.
%
%   MARKERS SHOW 514. The 10 sensitivity-flagged trials (d_final_cm < 1 cm |
%   v0_cm_s < 55 cm/s) are omitted FROM THE PLOTTED MARKERS ONLY, for visual
%   clarity. They remain in every fit. There is no separate sensitivity fit
%   and no separate marker style for them: a reader must not be able to
%   mistake the display selection for an analysis selection, because it is
%   not one.
%
%   So: a curve drawn through 514 visible points was fitted to 524. The 10
%   invisible trials pull that curve exactly as much as any other trial does.
%
%   NO ZERO-DROP CONTENT. These figures carry no v0 = 0 points, no measured
%   d_s markers, and read no d0_static export. The quasi-static d_s result is
%   reported in the text, not on these axes.
%
% ═══════════════════════════════════════════════════════════════════════════
%
%   USAGE
%       R = fig_scaling
%       R = fig_scaling('Root', 'D:\ME_GRANULAB\JerboaImpact')
%       R = fig_scaling('Save', false)
%
%   OPTIONS (name-value)
%       'Root'      results root  (default D:\ME_GRANULAB\JerboaImpact)
%       'Trials'    master trials CSV (default the approved 20260822_215312)
%       'OutDir'    output folder (default <Root>/03_RESULTS/_figures)
%       'Save'      write PDF + PNG                        (default true)
%       'Show'      display the figures                    (default true)
%
%   OUTPUTS (in OutDir)
%       fig_scaling.{pdf,png}              pooled 2/3 law
%       fig_scaling_geometry.{pdf,png}     (a) M1 per geometry (b) hull-normalised
%       fig_scaling_diagnostics.{pdf,png}  (a) residuals (b) RSS(n) profile
%
%   Base MATLAB only.

% ── options ──────────────────────────────────────────────────────────────
opt.Root   = 'D:\ME_GRANULAB\JerboaImpact';
opt.Trials = '';
opt.OutDir = '';
opt.Save   = true;
opt.Show   = true;
optNames = fieldnames(opt);
for i = 1:2:numel(varargin)
    j = find(strcmpi(optNames, varargin{i}), 1);
    if isempty(j)
        error('fig_scaling:unknownOption', ...
              'Unknown option "%s". Valid: %s', string(varargin{i}), strjoin(optNames', ', '));
    end
    opt.(optNames{j}) = varargin{i+1};
end
if isempty(opt.Trials)
    opt.Trials = fullfile(opt.Root, '03_RESULTS', '_exports', ...
                          'master_trials_20260822_215312.csv');
end
if isempty(opt.OutDir)
    opt.OutDir = fullfile(opt.Root, '03_RESULTS', '_figures');
end

% ── approved values (assertion targets; NOT fitted here) ─────────────────
APPROVED = struct( ...
    'd0',   0.044,  ...                     % pooled fixed-2/3 intercept, cm
    'a',    0.0711, ...                     % pooled fixed-2/3 coefficient
    'aG',   [0.0718 0.0703 0.0720], ...     % M1 slopes, order MDL below
    'nMin', 0.606,  ...                     % profile RSS argmin
    'nCI',  [0.20 0.99], ...                % cluster CI, approved; NOT re-bootstrapped
    'N',    524);
TOL3SF   = 5e-4;
TOL_NMIN = 0.005;

% Geometry identity. Colours are MATLAB's default order exactly as gscatter
% assigned it in step3_fits_working.png, where the groups sort as
% Default/Tight/Wide -> colours 1/2/3. Verified against that figure's legend.
% Marker shape is the one deliberate addition.
MDL  = ["Tight" "Default" "Wide"];
MARK = ["o" "s" "^"];
COL  = [0.8500 0.3250 0.0980;      % Tight   -> default colour 2 (orange/red)
        0      0.4470 0.7410;      % Default -> default colour 1 (blue)
        0.9290 0.6940 0.1250];     % Wide    -> default colour 3 (yellow)
MSZ  = 4;                          % matches the working figure's dot density
A_HULL = [2.607 3.495 4.052];      % cm^2, T/D/W; Fig 2(b) normalisation only

fprintf('\n=== fig_scaling ===\n');
fprintf('  Root   : %s\n', opt.Root);
fprintf('  Trials : %s\n', opt.Trials);
fprintf('  OutDir : %s\n', opt.OutDir);

% ── 1. trials and the approved selection ─────────────────────────────────
if ~isfile(opt.Trials)
    error('fig_scaling:noTrials', 'Master trials export not found:\n  %s', opt.Trials);
end
T = readtable(opt.Trials);
K = T(logical(T.keep_reviewed) & ~logical(T.isZeroDrop), :);
K.model = string(K.model);
N = height(K);
if N ~= APPROVED.N
    error('fig_scaling:wrongN', ...
        ['Selection gives N = %d, approved is %d. keep_reviewed & ~isZeroDrop is ' ...
         'the approved selection and this figure may not redefine it.'], N, APPROVED.N);
end
x   = K.v0_cm_s;
y   = K.d_final_cm;
x23 = x.^(2/3);

% Display-only subset. See CAPTION NOTES: this NEVER touches a fit.
flag  = K.d_final_cm < 1 | K.v0_cm_s < 55;
disp_ = ~flag;
nDisp = sum(disp_);

fprintf('\n--- selection ---\n');
fprintf('  fitted : N = %d  (keep_reviewed & ~isZeroDrop)\n', N);
fprintf('  shown  : n = %d  markers\n', nDisp);
fprintf('  Fits use N = %d; markers show %d (%d flagged trials omitted from display only).\n', ...
        N, nDisp, sum(flag));

% ── 2. pooled fixed-2/3 fit, on ALL 524 ──────────────────────────────────
cP = [ones(N,1) x23] \ y;
local_assert('pooled d0', cP(1), APPROVED.d0, TOL3SF);
local_assert('pooled a',  cP(2), APPROVED.a,  TOL3SF);

% ── 3. M1 per-geometry, shared d0, on ALL 524 ────────────────────────────
Gm = double(K.model == MDL);                 % N x 3, column order = MDL
c1 = [ones(N,1) Gm.*x23] \ y;                % [d0 aT aD aW]
for m = 1:3
    local_assert(sprintf('M1 a_%s', MDL(m)), c1(1+m), APPROVED.aG(m), TOL3SF);
end

% ── 4. profile RSS(n), recomputed exactly as step3_fits ──────────────────
ngrid = 0.2:0.01:1.4;
rssn  = arrayfun(@(nn) local_rss(x, y, nn), ngrid);
[~, iMin] = min(rssn);
nMin = ngrid(iMin);
local_assert('profile argmin n', nMin, APPROVED.nMin, TOL_NMIN);

fprintf('\n--- fits (recomputed on N = %d, asserted against approved) ---\n', N);
fprintf('  pooled fixed-2/3 : d0 = %.4f cm   a = %.4f   (approved %.3f / %.4f)\n', ...
        cP(1), cP(2), APPROVED.d0, APPROVED.a);
fprintf('  M1 shared d0     : d0 = %.4f cm\n', c1(1));
for m = 1:3
    fprintf('    a_%-8s     = %.4f         (approved %.4f)\n', MDL(m), c1(1+m), APPROVED.aG(m));
end
fprintf('  profile RSS(n)   : argmin n = %.3f (approved %.3f); CI [%.2f %.2f] indicated as approved\n', ...
        nMin, APPROVED.nMin, APPROVED.nCI);
fprintf('  RSS flatness     : RSS(2/3)/RSS_min = %.4f, RSS(1)/RSS_min = %.4f\n', ...
        local_rss(x,y,2/3)/rssn(iMin), local_rss(x,y,1)/rssn(iMin));

fprintf('\n--- geometry identity (fixed in every panel of every figure) ---\n');
fprintf('  %-8s %-6s %s\n', 'model', 'marker', 'colour (MATLAB default order)');
cname = ["default colour 2 (orange/red)" "default colour 1 (blue)" "default colour 3 (yellow)"];
for m = 1:3
    fprintf('  %-8s %-6s [%.4f %.4f %.4f]  %s\n', MDL(m), MARK(m), COL(m,:), cname(m));
end
fprintf('  no zero-drop content: no v0 = 0 points, no d_s markers, no d0_static read\n');

% ── 5. figures ───────────────────────────────────────────────────────────
S = struct('MDL',MDL,'MARK',MARK,'COL',COL,'MSZ',MSZ,'A_HULL',A_HULL, ...
           'x',x,'y',y,'model',K.model,'disp',disp_, ...
           'cP',cP,'c1',c1,'ngrid',ngrid,'rssn',rssn,'nMin',nMin, ...
           'nCI',APPROVED.nCI,'Show',opt.Show);

f1 = local_fig1(S);
f2 = local_fig2(S);
f3 = local_fig3(S);

written = strings(0,1);
if opt.Save
    written = [written; local_export(f1, opt.OutDir, 'fig_scaling')];
    written = [written; local_export(f2, opt.OutDir, 'fig_scaling_geometry')];
    written = [written; local_export(f3, opt.OutDir, 'fig_scaling_diagnostics')];
    fprintf('\n--- written ---\n');
    fprintf('  %s\n', written);
end
if ~opt.Show
    close(f1); close(f2); close(f3);
    f1 = gobjects(1); f2 = gobjects(1); f3 = gobjects(1);
end

fprintf('\n');
R = struct();
R.N          = N;
R.nDisplayed = nDisp;
R.pooled     = cP;
R.M1         = c1;
R.nMin       = nMin;
R.nCI        = APPROVED.nCI;
R.ngrid      = ngrid;
R.rssn       = rssn;
R.written    = written;
R.figures    = [f1 f2 f3];
end

% ═════════════════════════════════════════════════════════════════════════
%  FIGURE 1 -- pooled 2/3 law, styled as the step3 working left panel
% ═════════════════════════════════════════════════════════════════════════
function fig = local_fig1(S)
    fig = figure('Color','w','Units','inches','Position',[1 1 3.375 2.85], ...
                 'Visible', local_tern(S.Show,'on','off'));
    ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    h = gobjects(0); lbl = strings(0,1);
    for m = 1:3
        r = S.disp & S.model == S.MDL(m);
        h(end+1) = local_pts(ax, S.x(r), S.y(r), S, m); %#ok<AGROW>
        lbl(end+1,1) = S.MDL(m); %#ok<AGROW>
    end

    % Pooled curve over the span of the PLOTTED data only: no extrapolation,
    % and in particular no extension toward v0 = 0, which these data do not
    % reach and this figure does not claim.
    vv = linspace(min(S.x(S.disp)), max(S.x(S.disp)), 200);
    h(end+1) = plot(ax, vv, S.cP(1) + S.cP(2)*vv.^(2/3), 'k-', 'LineWidth', 1.0);
    lbl(end+1,1) = "fit";

    xlabel(ax, 'v_0 (cm/s)');
    ylabel(ax, 'd (cm)');
    legend(ax, h, cellstr(lbl), 'Location', 'southeast');
end

% ═════════════════════════════════════════════════════════════════════════
%  FIGURE 2 -- per-geometry M1, raw and hull-normalised
% ═════════════════════════════════════════════════════════════════════════
function fig = local_fig2(S)
%   Same 1x2 orientation as the step4 working figure, restyled to the step3
%   look and the geometry colours. No per-height overlays: the individual
%   trials are the evidence.
    fig = figure('Color','w','Units','inches','Position',[1 1 6.75 2.85], ...
                 'Visible', local_tern(S.Show,'on','off'));
    tl = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');
    vv = linspace(min(S.x(S.disp)), max(S.x(S.disp)), 200);
    tags = ["(a)" "(b)"];
    for p = 1:2
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for m = 1:3
            s = local_tern(p == 1, 1, S.A_HULL(m)^(1/3));
            r = S.disp & S.model == S.MDL(m);
            local_pts(ax, S.x(r), S.y(r)./s, S, m);
        end
        hh = gobjects(1,3);
        for m = 1:3
            s = local_tern(p == 1, 1, S.A_HULL(m)^(1/3));
            hh(m) = plot(ax, vv, (S.c1(1) + S.c1(1+m)*vv.^(2/3))./s, '-', ...
                'Color', S.COL(m,:), 'LineWidth', 1.4);
        end
        xlabel(ax, 'v_0 (cm/s)');
        if p == 1
            ylabel(ax, 'd (cm)');
            legend(ax, hh, cellstr(S.MDL), 'Location', 'southeast');
        else
            ylabel(ax, 'd / A_{hull}^{1/3}');
        end
        local_panel_tag(ax, tags(p));
    end
end

% ═════════════════════════════════════════════════════════════════════════
%  FIGURE 3 -- diagnostics
% ═════════════════════════════════════════════════════════════════════════
function fig = local_fig3(S)
    fig = figure('Color','w','Units','inches','Position',[1 1 6.75 2.85], ...
                 'Visible', local_tern(S.Show,'on','off'));
    tl = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

    % ── (a) residuals, as the step3 working middle panel ─────────────────
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    res = S.y - (S.cP(1) + S.cP(2)*S.x.^(2/3));
    yline(ax, 0, 'k-');
    hh = gobjects(1,3);
    for m = 1:3
        r = S.disp & S.model == S.MDL(m);
        hh(m) = local_pts(ax, S.x(r), res(r), S, m);
    end
    xlabel(ax, 'v_0 (cm/s)');
    ylabel(ax, 'residual (cm)');
    legend(ax, hh, cellstr(S.MDL), 'Location', 'southwest');
    local_panel_tag(ax, "(a)");

    % ── (b) the RSS profile, PRESERVED from the working right panel ──────
    %   Same curve appearance (default line, default colour), same reference
    %   line styles, same autoscaled ranges. The only addition is the CI band;
    %   the only removals are the working title and its inline CI text, which
    %   the caption carries instead.
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, S.ngrid, S.rssn/min(S.rssn), '-');
    yl = ylim(ax);                    % the curve alone sets the y range, as before
    pb = patch(ax, [S.nCI(1) S.nCI(2) S.nCI(2) S.nCI(1)], [yl(1) yl(1) yl(2) yl(2)], ...
               [0.88 0.88 0.88], 'EdgeColor', 'none');
    uistack(pb, 'bottom');            % behind the curve and the grid lines
    xline(ax, 2/3, '--');
    xline(ax, S.nMin);
    xline(ax, 1, ':');
    % Pinned, not autoscaled: the band and the reference lines must not widen
    % the panel away from the working figure's ranges.
    xlim(ax, [S.ngrid(1) S.ngrid(end)]);
    ylim(ax, yl);
    xlabel(ax, 'n');
    ylabel(ax, 'RSS(n)/RSS_{min}');
    local_panel_tag(ax, "(b)");
end

% ═════════════════════════════════════════════════════════════════════════
%  helpers
% ═════════════════════════════════════════════════════════════════════════
function h = local_pts(ax, xx, yy, S, m)
%LOCAL_PTS  One geometry's trials, in the working figure's density: opaque
%   filled markers at the default colour, differing only in shape.
    h = plot(ax, xx, yy, 'LineStyle', 'none', ...
             'Marker', S.MARK(m), 'MarkerSize', S.MSZ, ...
             'MarkerFaceColor', S.COL(m,:), 'MarkerEdgeColor', S.COL(m,:));
end

function rss = local_rss(x, y, n)
%LOCAL_RSS  RSS of d = d0 + a*v0^n, the same least squares step3_fits.m uses,
%   so the profile here cannot drift from the analysis.
    X = [ones(size(x)) x.^n];
    rss = sum((y - X*(X\y)).^2);
end

function local_assert(what, got, want, tol)
%LOCAL_ASSERT  A recomputed coefficient must match its approved value.
    if ~isfinite(got) || abs(got - want) > tol
        error('fig_scaling:assertFailed', ...
            ['%s recomputed as %.6g, approved %.6g (tolerance %.1g). This figure ' ...
             'reproduces approved results and does not refit: a mismatch means the ' ...
             'inputs or the selection changed, and the figure must not be drawn ' ...
             'until that is understood.'], what, got, want, tol);
    end
end

function local_panel_tag(ax, s)
%LOCAL_PANEL_TAG  Panel letter inside the top-left corner. Inside on purpose:
%   a tiledlayout clips anything placed outside the axes, so an outside tag
%   silently disappears from the exported file.
    text(ax, 0.03, 0.97, s, 'Units','normalized', 'FontWeight','bold', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top');
end

function paths = local_export(fig, outDir, stem)
%LOCAL_EXPORT  Vector PDF for the manuscript, 300-dpi PNG alongside.
    if ~isfolder(outDir), mkdir(outDir); end
    pdfPath = fullfile(outDir, [stem '.pdf']);
    pngPath = fullfile(outDir, [stem '.png']);
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    exportgraphics(fig, pngPath, 'Resolution', 300);
    paths = [string(pdfPath); string(pngPath)];
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
