function R = fig_foot_schematic(varargin)
%FIG_FOOT_SCHEMATIC  Publication schematics of the three foot models, from CAD.
%
%   READ-ONLY with respect to the results tree: it reads the STLs in cad/ and
%   writes only figures and one CSV to OutDir.
%
%   Every schematic is RENDERED FROM THE GEOMETRY. Nothing is drawn by hand,
%   traced, or approximated, so a figure cannot drift from the model it claims
%   to show. If an STL is missing the script errors naming the path rather than
%   substituting anything.
%
%   USAGE
%       fig_foot_schematic
%       fig_foot_schematic('Root', 'D:\ME_GRANULAB\JerboaImpact')
%       R = fig_foot_schematic('Save', false, 'Labels', true);
%
%   OPTIONS (name-value)
%       'CadDir'  STL folder                (default <repo>/cad)
%       'CutY'    cut plane, mm             (default -50, the beam top)
%       'Root'    results root; when given, OutDir defaults to
%                 <Root>/03_RESULTS/_figures
%       'OutDir'  output folder             (default per 'Root', else ./figures)
%       'Save'    write PDF + PNG + CSV     (default true)
%       'Show'    display the figures       (default true)
%       'Labels'  add callouts to Figure A  (default false)
%       'ScaleBar' draw a 10 mm bar on Figure A (default false). The axes are
%                 on and ticked in mm, so a bar duplicates the scale and can
%                 disagree with it after a resize. Kept only for slides.
%
%   WHAT IS COMPUTED
%     A_bare   area of the true projected footprint at y <= CutY, onto XZ
%              (normal to the drop axis). Holes between the toes are real holes:
%              polyshape keeps them, so they are not counted as intruding area.
%     A_hull   area of the convex hull of the same projected points. The
%              envelope the foot sweeps, ignoring the gaps.
%
%   The three models differ ONLY in toe splay, so A_bare is identical for all
%   three and only A_hull grows. That is the point of the figure.
%
%   VALIDATION. Both areas are asserted against values computed independently
%   (trimesh/shapely, 2026-08-20) at 2% relative tolerance, on every run. A
%   mismatch aborts before any figure is written: a schematic that disagrees
%   with the reference geometry must not reach a paper.
%
%   Base MATLAB only -- stlread, triangulation, polyshape, convhull,
%   exportgraphics. No additional toolboxes.

% ── options ──────────────────────────────────────────────────────────────
repoRoot = fileparts(fileparts(mfilename('fullpath')));

opt.CadDir = fullfile(repoRoot, 'cad');
opt.CutY   = -50;
opt.Root   = '';
opt.OutDir = '';
opt.Save   = true;
opt.Show   = true;
opt.Labels   = false;
opt.ScaleBar = false;   % axes ticks carry the scale; kept for the odd talk slide
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

if isempty(opt.OutDir)
    if ~isempty(opt.Root)
        opt.OutDir = fullfile(opt.Root, '03_RESULTS', '_figures');
    else
        opt.OutDir = fullfile(pwd, 'figures');
    end
end

% Model order is deliberate: increasing toe splay, left to right.
MODELS = { ...
    'Tight',   'jerboa_foot_model_tighter_toe_spacing.stl' ; ...
    'Default', 'jerboa_foot_model_rectangularbeam.stl'     ; ...
    'Wide',    'jerboa_foot_model_wider_toe_spacing.stl'   };

% Independently computed reference areas, cm^2 (trimesh/shapely, 2026-08-20).
REF_BARE = [2.122; 2.122; 2.122];
REF_HULL = [2.607; 3.495; 4.052];
REF_CUTY = -50;
TOL_REL  = 0.02;

fprintf('\n=== fig_foot_schematic ===\n');
fprintf('  CadDir : %s\n', opt.CadDir);
fprintf('  CutY   : %g mm (beam top; everything below this enters the bed)\n', opt.CutY);
fprintf('  OutDir : %s\n', opt.OutDir);
fprintf('  Save   : %s\n', local_tern(opt.Save,'yes','no'));
fprintf('  Show   : %s\n', local_tern(opt.Show,'yes','no'));

% ── geometry ─────────────────────────────────────────────────────────────
n = size(MODELS,1);
name    = strings(n,1);
stlPath = strings(n,1);
Abare   = nan(n,1);
Ahull   = nan(n,1);
FP      = cell(n,1);      % footprint polyshape (clipped, XZ)
HULL    = cell(n,1);      % hull polygon [x z]
TRI     = cell(n,1);      % full triangulation, reused by Figures A and C

fprintf('\n--- geometry ---\n');
for i = 1:n
    name(i)    = string(MODELS{i,1});
    stlPath(i) = string(fullfile(opt.CadDir, MODELS{i,2}));

    if ~isfile(stlPath(i))
        error('fig_foot_schematic:missingSTL', ...
            ['Required STL not found:\n  %s\n\nThe schematics are rendered ' ...
             'from the CAD geometry and nothing is substituted for a missing ' ...
             'file. Export it from the OpenSCAD model into %s -- see ' ...
             'cad/README.md for the expected filenames.'], ...
            stlPath(i), opt.CadDir);
    end

    % Read once and keep: Figures A and C project the FULL mesh, and reading
    % the file again there could silently pick up a different revision from the
    % one validation just ran on.
    tri    = stlread(char(stlPath(i)));
    TRI{i} = tri;
    tris   = local_clip_below(tri, opt.CutY);
    if isempty(tris)
        error('fig_foot_schematic:nothingBelowCut', ...
            ['%s has no geometry at y <= %g. Check CutY against the STL ' ...
             'convention in cad/README.md (drop axis is -Y, toe tips near ' ...
             'y = -61.2).'], MODELS{i,2}, opt.CutY);
    end

    [FP{i}, Abare(i)]   = local_footprint(tris);
    [HULL{i}, Ahull(i)] = local_hull(tris);

    fprintf('  %-8s %5d clipped tri  A_bare = %.3f cm^2  A_hull = %.3f cm^2  (hull/bare %.2f)\n', ...
            name(i), size(tris,1), Abare(i), Ahull(i), Ahull(i)/Abare(i));
end

% ── validation ───────────────────────────────────────────────────────────
% Runs before anything is written. A schematic that disagrees with the
% independently computed geometry must not reach a figure file.
fprintf('\n--- validation (2%% relative, vs trimesh/shapely 2026-08-20) ---\n');
if abs(opt.CutY - REF_CUTY) > 1e-9
    fprintf(['  SKIPPED: reference areas are defined at CutY = %g, this run ' ...
             'used %g.\n           Areas below are still computed from the ' ...
             'geometry, but are unvalidated.\n'], REF_CUTY, opt.CutY);
else
    ok = true;
    for i = 1:n
        [okB, msgB] = local_check(name(i), 'A_bare', Abare(i), REF_BARE(i), TOL_REL);
        [okH, msgH] = local_check(name(i), 'A_hull', Ahull(i), REF_HULL(i), TOL_REL);
        fprintf('%s%s', msgB, msgH);
        ok = ok && okB && okH;
    end
    if ~ok
        error('fig_foot_schematic:validationFailed', ...
            ['Computed areas disagree with the reference values by more than ' ...
             '%g%%. No figure was written.\nEither the CAD was revised (update ' ...
             'the reference values in this script AND cad/README.md, saying ' ...
             'why), or the geometry pipeline has changed. Do not publish a ' ...
             'figure until the two agree.'], 100*TOL_REL);
    end
    fprintf('  all within tolerance\n');
end

% ── figures ──────────────────────────────────────────────────────────────
if opt.Save && ~isfolder(opt.OutDir), mkdir(opt.OutDir); end

iDefault = find(name == "Default", 1);
figA = local_figure_assembly(TRI{iDefault}, opt);
figB = local_figure_footprints(name, FP, HULL, Abare, Ahull, opt);
figC = local_figure_topview(name, TRI, HULL, opt);

written = strings(0,1);
if opt.Save
    written = [written; local_export(figA, opt.OutDir, 'fig_foot_assembly')];
    written = [written; local_export(figB, opt.OutDir, 'fig_foot_footprints')];
    written = [written; local_export(figC, opt.OutDir, 'fig_foot_topview')];
    csvPath = local_write_csv(opt, name, stlPath, Abare, Ahull);
    written = [written; string(csvPath)];
    fprintf('\n--- written ---\n');
    fprintf('  %s\n', written);
end

if ~opt.Show
    close(figA); close(figB); close(figC);
    figA = gobjects(1); figB = gobjects(1); figC = gobjects(1);
end

fprintf('\n');
% Built field by field: struct() with array values would return a 1-by-3
% struct array rather than one struct holding the three models.
R = struct();
R.model          = name;
R.stl            = stlPath;
R.A_bare_cm2     = Abare;
R.A_hull_cm2     = Ahull;
R.hull_over_bare = Ahull ./ Abare;
R.CutY           = opt.CutY;
R.footprint      = FP;
R.hull           = HULL;
R.written        = written;
R.figAssembly    = figA;
R.figFootprints  = figB;
R.figTopView     = figC;
end

% ═════════════════════════════════════════════════════════════════════════
%  GEOMETRY CORE
% ═════════════════════════════════════════════════════════════════════════
function tris = local_clip_below(tri, ycut)
%LOCAL_CLIP_BELOW  Clip every triangle against the half-space y <= ycut.
%
%   Sutherland-Hodgman, per triangle. A triangle crossing the plane yields a
%   3- or 4-gon; 4-gons are fan-triangulated. Triangles are NOT simply dropped
%   when they cross: doing that would bite chunks out of the footprint exactly
%   at the cut, which is where the area is being measured.
%
%   Returns an N-by-9 matrix, one triangle per row:
%       [x1 y1 z1  x2 y2 z2  x3 y3 z3]
%   Each output triangle carries its own vertices; duplicates across triangles
%   are harmless for both the polyshape union and the convex hull.

    P = tri.Points;
    F = tri.ConnectivityList;
    nF = size(F,1);

    % Worst case a crossing triangle becomes two, so preallocate 2*nF.
    tris = zeros(2*nF, 9);
    k = 0;

    for f = 1:nF
        V = P(F(f,:), :);              % 3x3, rows are vertices
        poly = local_clip_poly(V, ycut);
        m = size(poly,1);
        if m < 3, continue; end
        % Fan-triangulate: 3-gon gives one triangle, 4-gon gives two.
        for t = 2:(m-1)
            k = k + 1;
            tris(k,:) = [poly(1,:), poly(t,:), poly(t+1,:)];
        end
    end
    tris = tris(1:k, :);
end

function out = local_clip_poly(V, ycut)
%LOCAL_CLIP_POLY  Sutherland-Hodgman clip of one convex polygon to y <= ycut.
%   V is m-by-3. Walks the edges, keeping inside vertices and inserting the
%   crossing point on every edge that changes side.
    m = size(V,1);
    out = zeros(0,3);
    for i = 1:m
        A = V(i,:);
        B = V(mod(i, m) + 1, :);
        aIn = A(2) <= ycut;
        bIn = B(2) <= ycut;
        if aIn
            out(end+1,:) = A; %#ok<AGROW>
        end
        if aIn ~= bIn
            % Edge crosses the plane; t is exact because the plane is axis-aligned.
            dy = B(2) - A(2);
            if abs(dy) > 0             % guaranteed non-zero when sides differ
                t = (ycut - A(2)) / dy;
                out(end+1,:) = A + t*(B - A); %#ok<AGROW>
            end
        end
    end
end

function [pg, areaCm2] = local_footprint(tris)
%LOCAL_FOOTPRINT  Union of the clipped triangles projected onto XZ.
%
%   XZ is normal to the drop axis, so this is the true intruding footprint.
%   Holes between the toes are represented as holes by polyshape and are
%   therefore NOT counted as intruding area -- which is the whole distinction
%   between this and the convex hull.

    % polyshape repairs degenerate input loudly; that is expected here, since a
    % triangle mesh has many near-collinear projections. Silence locally only.
    ws = warning;                     % full current state, restored on exit
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    warning('off', 'MATLAB:polyshape:boundary3Points');
    restore = onCleanup(@() warning(ws)); %#ok<NASGU>

    nT = size(tris,1);
    parts = cell(nT,1);
    kept  = 0;
    for i = 1:nT
        x = tris(i, [1 4 7]);
        z = tris(i, [3 6 9]);
        if polyarea(x, z) < 1e-9       % edge-on triangle: no projected area
            continue
        end
        kept = kept + 1;
        parts{kept} = polyshape(x, z, 'Simplify', true);
    end
    if kept == 0
        error('fig_foot_schematic:emptyFootprint', ...
            ['Every clipped triangle projected to zero area on XZ. The cut ' ...
             'plane or the STL axis convention is wrong -- see cad/README.md.']);
    end
    parts = parts(1:kept);
    pg = union([parts{:}]);

    areaCm2 = area(pg) / 100;          % mm^2 -> cm^2
end

function [hullXZ, areaCm2] = local_hull(tris)
%LOCAL_HULL  Convex hull of ALL clipped vertices projected onto XZ.
%   The envelope swept by the foot, ignoring the gaps between the toes.
    X = tris(:, [1 4 7]); X = X(:);
    Z = tris(:, [3 6 9]); Z = Z(:);
    k = convhull(X, Z);
    hullXZ  = [X(k), Z(k)];
    areaCm2 = polyarea(X(k), Z(k)) / 100;
end

function [ok, msg] = local_check(model, what, got, ref, tolRel)
    rel = abs(got - ref) / abs(ref);
    ok  = rel <= tolRel;
    msg = sprintf('  %-8s %-7s computed %.3f  reference %.3f  (%.2f%%)  %s\n', ...
                  model, what, got, ref, 100*rel, local_tern(ok,'ok','FAIL'));
end

% ═════════════════════════════════════════════════════════════════════════
%  FIGURE A -- assembly side view
% ═════════════════════════════════════════════════════════════════════════
function pg = local_project_union(tri, iAx, jAx)
%LOCAL_PROJECT_UNION  Silhouette of a FULL mesh projected onto two axes.
%
%   Shared by Figure A (project along Z onto XY) and Figure C (project along Y
%   onto XZ), so the two views are built by the same code and cannot drift into
%   different renderings of the same part.
%
%   iAx/jAx are column indices into tri.Points: 1 = x, 2 = y, 3 = z.
%   Returns the polyshape union; holes survive as holes.

    ws = warning;                     % full current state, restored on exit
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    warning('off', 'MATLAB:polyshape:boundary3Points');
    restore = onCleanup(@() warning(ws)); %#ok<NASGU>

    P = tri.Points;  F = tri.ConnectivityList;
    parts = cell(size(F,1),1);
    kept = 0;
    for f = 1:size(F,1)
        u = P(F(f,:), iAx).';
        v = P(F(f,:), jAx).';
        if polyarea(u, v) < 1e-9, continue; end   % edge-on: no projected area
        kept = kept + 1;
        parts{kept} = polyshape(u, v, 'Simplify', true);
    end
    if kept == 0
        error('fig_foot_schematic:emptyProjection', ...
            'Every triangle projected to zero area on axes %d/%d.', iAx, jAx);
    end
    pg = union([parts{1:kept}]);
end

function local_style_axes(ax, xlab, ylab, ttl, fs)
%LOCAL_STYLE_AXES  The one axes convention every panel in this file uses.
%   Engineering orthographic: visible box, equal aspect, mm ticks, labelled
%   axes. Kept in one place so the figures cannot drift apart stylistically.
    box(ax, 'on');
    daspect(ax, [1 1 1]);     % equal aspect that leaves explicit limits alone
    set(ax, 'FontSize', fs, 'LineWidth', 0.5, 'Layer', 'top', ...
            'XGrid','off', 'YGrid','off');
    xlabel(ax, xlab, 'FontSize', fs);
    ylabel(ax, ylab, 'FontSize', fs);
    if ~isempty(ttl), title(ax, ttl, 'FontSize', fs+1); end
end

% ═════════════════════════════════════════════════════════════════════════
%  FIGURE A -- assembly side view
% ═════════════════════════════════════════════════════════════════════════
function fig = local_figure_assembly(tri, opt)
%LOCAL_FIGURE_ASSEMBLY  Default model, FULL mesh, projected along Z onto XY.
%   Silhouette only -- the polyshape union filled grey with a thin boundary, no
%   wireframe. Mesh lines would read as structure that is not there.

    pg = local_project_union(tri, 1, 2);      % x, y

    % Data is roughly 100 x 80 mm, so a ~6 x 5 in canvas holds it at equal
    % aspect with no dead whitespace.
    fig = figure('Color','w','Units','inches','Position',[1 1 6.0 5.0], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    ax = axes(fig); hold(ax,'on');

    plot(pg, 'Parent', ax, 'FaceColor', [0.85 0.85 0.85], 'FaceAlpha', 1, ...
             'EdgeColor', 'k', 'LineWidth', 0.8);

    % Limits from the geometry's own bounding box plus a fixed 5 mm margin,
    % rather than autoscale: reproducible, and no dead whitespace.
    [bx, by] = boundingbox(pg);
    MARGIN = 5;
    xl = [bx(1)-MARGIN, bx(2)+MARGIN];
    yl = [by(1)-MARGIN, by(2)+MARGIN];
    xlim(ax, xl);  ylim(ax, yl);

    % ── cut plane ────────────────────────────────────────────────────────
    plot(ax, xl, [opt.CutY opt.CutY], 'k--', 'LineWidth', 1.0);
    % Label at the LEFT end, above the line. At y = CutY the part lies to the
    % right, so the left end is clear of the silhouette.
    text(ax, xl(1) + 0.02*diff(xl), opt.CutY, ' area cut (beam top)', ...
         'FontSize', 8, 'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'bottom');

    % ── drop direction ───────────────────────────────────────────────────
    % Short arrow inside the axes, pointing -Y. Its column is CHOSEN by testing
    % candidates against the silhouette, so it lands in genuinely free space
    % rather than at a coordinate that happened to be clear for one CAD
    % revision.
    ARROW_LEN = 15;                                   % mm
    yTop   = yl(2) - 0.06*diff(yl);
    xArrow = local_clear_column(pg, xl, yTop, ARROW_LEN);
    quiver(ax, xArrow, yTop, 0, -ARROW_LEN, 'Color', 'k', 'LineWidth', 1.1, ...
           'MaxHeadSize', 0.6, 'AutoScale', 'off');
    text(ax, xArrow + 0.015*diff(xl), yTop - ARROW_LEN/2, 'drop', ...
         'FontSize', 8, 'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'middle');

    % ── optional scale bar ───────────────────────────────────────────────
    if opt.ScaleBar
        % Off by default: the axes are ticked in mm, so a bar restates the
        % scale and can contradict it after a resize.
        sbX = xl(2) - 0.05*diff(xl) - 10;
        sbY = yl(1) + 0.06*diff(yl);
        plot(ax, [sbX sbX+10], [sbY sbY], 'k-', 'LineWidth', 2.0);
        text(ax, sbX+5, sbY, '10 mm', 'FontSize', 8, ...
             'HorizontalAlignment','center', 'VerticalAlignment','top');
    end

    if opt.Labels
        % Positions derived from the model's own y-extent, not typed in, so
        % they follow the geometry if the CAD changes.
        yspan = diff(yl);
        xr2   = xl(2) - 0.02*diff(xl);
        lab = { 'mount',   yl(2) - 0.08*yspan ; ...
                'linkage', yl(2) - 0.30*yspan ; ...
                'beam',    yl(2) - 0.62*yspan ; ...
                'toes',    yl(1) + 0.06*yspan };
        for i = 1:size(lab,1)
            text(ax, xr2, lab{i,2}, lab{i,1}, 'FontSize', 8, ...
                 'HorizontalAlignment','right', 'VerticalAlignment','middle');
        end
    end

    local_style_axes(ax, 'x (mm)', 'y (mm)', ...
                     'Foot assembly (Default), viewed along z', 8);
end

function xBest = local_clear_column(pg, xl, yTop, arrowLen)
%LOCAL_CLEAR_COLUMN  An x where a vertical arrow of arrowLen misses the part.
%
%   Scans candidate columns from the left and returns the first whose whole
%   span is outside the silhouette. Falls back to the far left if the part
%   spans everything, which at worst overlaps rather than erroring: a figure
%   should still be produced.
    ys = linspace(yTop - arrowLen, yTop, 12);
    for f = 0.04 : 0.02 : 0.50
        x = xl(1) + f*diff(xl);
        if ~any(isinterior(pg, repmat(x, size(ys)), ys))
            xBest = x;
            return
        end
    end
    xBest = xl(1) + 0.04*diff(xl);
end

% ═════════════════════════════════════════════════════════════════════════
%  FIGURE B -- intruding footprints
% ═════════════════════════════════════════════════════════════════════════
function fig = local_figure_footprints(name, FP, HULL, Abare, Ahull, opt)
%LOCAL_FIGURE_FOOTPRINTS  1x3, increasing toe splay, shared limits.
%   Shared x/z limits are what make the panels comparable: the growth of the
%   hull across models is the message, and per-panel autoscaling would hide it.

    n = numel(name);
    fig = figure('Color','w','Units','inches','Position',[1 1 7.0 2.6], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl = tiledlayout(fig, 1, n, 'Padding','compact', 'TileSpacing','compact');

    % Common limits over every footprint AND hull, with a small margin.
    xs = []; zs = [];
    for i = 1:n
        [xi, zi] = boundingbox(FP{i});
        xs = [xs, xi]; zs = [zs, zi]; %#ok<AGROW>
        xs = [xs, min(HULL{i}(:,1)), max(HULL{i}(:,1))]; %#ok<AGROW>
        zs = [zs, min(HULL{i}(:,2)), max(HULL{i}(:,2))]; %#ok<AGROW>
    end
    MARGIN = 2;                      % mm, fixed: comparable panels, no dead space
    XL = [min(xs)-MARGIN, max(xs)+MARGIN];
    ZL = [min(zs)-MARGIN, max(zs)+MARGIN];

    for i = 1:n
        ax = nexttile(tl); hold(ax,'on'); box(ax,'on');

        % Bare footprint. polyshape renders holes as background, so the gaps
        % between the toes read as gaps rather than as filled area.
        plot(FP{i}, 'Parent', ax, 'FaceColor', [0.75 0.75 0.75], ...
             'FaceAlpha', 1, 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);

        % Convex hull
        plot(ax, HULL{i}(:,1), HULL{i}(:,2), 'r--', 'LineWidth', 1.2);

        % Shared limits, set before the style helper. Ticks and labels go on
        % EVERY panel: each is a standalone orthographic view, and a reader
        % should not have to look left to find the units.
        xlim(ax, XL); ylim(ax, ZL);
        local_style_axes(ax, 'x (mm)', 'z (mm)', '', 8);

        % Title carries the COMPUTED values, never typed-in ones.
        title(ax, sprintf('%s: A_{bare} = %.2f, A_{hull} = %.2f cm^2', ...
                          name(i), Abare(i), Ahull(i)), ...
              'FontSize', 8, 'FontWeight','normal');
    end

    title(tl, sprintf('Intruding footprint at y = %g mm (grey) and convex hull (red)', ...
                      opt.CutY), 'FontSize', 9);
end

% ═════════════════════════════════════════════════════════════════════════
%  FIGURE C -- plan views along the drop axis
% ═════════════════════════════════════════════════════════════════════════
function fig = local_figure_topview(name, TRI, HULL, opt)
%LOCAL_FIGURE_TOPVIEW  Full assembly seen along the drop axis, per model.
%
%   The companion to Figure A: A is the full assembly from the side, this is
%   the full assembly from the direction it travels. FULL mesh, no cut --
%   Figure B is the cut one.
%
%   The intruding footprint's convex hull is overlaid in the SAME red dashed
%   style Figure B uses, so the reader can see which part of the plan view
%   actually enters the bed. The hull polygons are the ones already computed
%   for Figure B and passed in; recomputing them here, even by the same method,
%   would risk the two figures disagreeing about the same quantity.

    n = numel(name);
    fig = figure('Color','w','Units','inches','Position',[1 1 7.0 2.8], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl = tiledlayout(fig, 1, n, 'Padding','compact', 'TileSpacing','compact');

    % Project every model first, so the shared limits cover all three. Same
    % helper Figure A uses, with the other axis pair: along Y onto XZ.
    PG = cell(n,1);
    xs = []; zs = [];
    for i = 1:n
        PG{i} = local_project_union(TRI{i}, 1, 3);     % x, z
        [xi, zi] = boundingbox(PG{i});
        xs = [xs, xi]; zs = [zs, zi]; %#ok<AGROW>
        xs = [xs, min(HULL{i}(:,1)), max(HULL{i}(:,1))]; %#ok<AGROW>
        zs = [zs, min(HULL{i}(:,2)), max(HULL{i}(:,2))]; %#ok<AGROW>
    end
    MARGIN = 2;                      % mm, fixed, as in Figure B
    XL = [min(xs)-MARGIN, max(xs)+MARGIN];
    ZL = [min(zs)-MARGIN, max(zs)+MARGIN];

    for i = 1:n
        ax = nexttile(tl); hold(ax,'on');

        % Same silhouette treatment as every other panel in this file: grey
        % fill, thin black boundary, holes left white, no wireframe.
        plot(PG{i}, 'Parent', ax, 'FaceColor', [0.85 0.85 0.85], ...
             'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 0.8);

        % The intruding region's hull, in Figure B's red dashed style.
        plot(ax, HULL{i}(:,1), HULL{i}(:,2), 'r--', 'LineWidth', 1.2);

        xlim(ax, XL); ylim(ax, ZL);
        local_style_axes(ax, 'x (mm)', 'z (mm)', '', 8);
        title(ax, name(i), 'FontSize', 8, 'FontWeight', 'normal');
    end

    title(tl, sprintf(['Full assembly viewed along the drop axis (grey), ' ...
                       'with the intruding hull at y = %g mm (red)'], opt.CutY), ...
          'FontSize', 9);
end

% ═════════════════════════════════════════════════════════════════════════
function paths = local_export(fig, outDir, stem)
%LOCAL_EXPORT  Vector PDF for the paper, 600-dpi PNG for everything else.
    pdfPath = fullfile(outDir, [stem '.pdf']);
    pngPath = fullfile(outDir, [stem '.png']);
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    exportgraphics(fig, pngPath, 'Resolution', 600);
    paths = [string(pdfPath); string(pngPath)];
end

function csvPath = local_write_csv(opt, name, stlPath, Abare, Ahull)
%LOCAL_WRITE_CSV  The computed areas, with enough provenance to be re-derived.
    csvPath = fullfile(opt.OutDir, 'foot_areas_computed.csv');
    fid = fopen(csvPath, 'w');
    if fid < 0
        error('fig_foot_schematic:csvFailed', 'Could not write %s', csvPath);
    end
    fprintf(fid, '# computed by scripts/fig_foot_schematic.m on %s\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# CutY = %g mm\n', opt.CutY);
    fprintf(fid, '# areas are COMPUTED from the STL geometry, not transcribed\n');
    for i = 1:numel(name)
        fprintf(fid, '# %-8s %s\n', name(i), stlPath(i));
    end
    fprintf(fid, 'model,A_bare_cm2,A_hull_cm2,hull_over_bare,CutY_mm\n');
    for i = 1:numel(name)
        fprintf(fid, '%s,%.4f,%.4f,%.4f,%g\n', ...
                name(i), Abare(i), Ahull(i), Ahull(i)/Abare(i), opt.CutY);
    end
    fclose(fid);
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
