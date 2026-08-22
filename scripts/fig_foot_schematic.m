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
%       'DepthSweep'  also compute A_bare(z), A_hull(z), a_local(z) over a
%                 depth grid, draw Figure C and write foot_area_vs_depth.csv
%                 (default true)
%       'SweepZ'  depth grid, mm below the toe tip (default [0:0.25:12 12.5:0.5:40])
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
%   DEPTH SWEEP (Figure C). The single-cut numbers above describe the FOOT
%   ONLY: everything at y <= CutY, i.e. the distal toe segments (y = -61.19
%   ... -53.20) and the base block + proximal toe bars (y = -53.20 ... -50.17).
%   The inclined bar ("linkage"), the marker post and the mount all start at
%   the y = -50 junction and are NOT in those numbers -- but they enter the bed
%   once the toe tip is deeper than z = 11.19 mm, which most trials exceed.
%   Figure C therefore reports, as a function of toe-tip depth z:
%     A_bare(z), A_hull(z)  cumulative projected area of EVERYTHING below the
%                           surface (same clip + projection as above)
%     a_local(z)            true cross-section of the model at the surface
%                           plane, from the signed projected area of the
%                           clipped open surface (divergence theorem; exact
%                           for a watertight, consistently wound STL)
%   Geometry-derived landmarks (toe tip, bar bottom, cut, post top, mount
%   bottom) are asserted against the documented values and drawn on Figures
%   A and C, so the labels cannot drift from the CAD.
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
opt.Labels = false;
opt.DepthSweep = true;
opt.SweepZ = [0:0.25:12, 12.5:0.5:40];   % mm below the toe tip
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

% Depth-sweep reference (trimesh/shapely, 2026-08-22), cm^2, at toe-tip depth
% z (mm). A_bare and a_local are identical across models to < 0.05 %; A_hull
% is per model [Tight Default Wide]. Asserted at 2 % when the sweep runs.
REF_SWEEP_Z    = [ 6     10     20     29   ];
REF_SWEEP_BARE = [0.594  2.014  2.957  3.729];
REF_SWEEP_LOC  = [0.120  1.392  0.388  0.388];
REF_SWEEP_HULL = [0.893  2.501  3.491  4.285 ;   % Tight
                  1.452  3.390  4.472  5.398 ;   % Default
                  1.803  3.946  5.150  6.195 ];  % Wide

% Documented landmarks (cad/README.md), mm. Re-derived from each STL by
% local_landmarks and asserted to 0.05 mm.
REF_LANDMARKS = struct('toeTip', -61.192, 'barBottom', -53.199, ...
                       'footTop', -50.173, 'postTop', -21.423, 'mountBottom', -2.5);

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
FP      = cell(n,1);      % footprint polyshape
HULL    = cell(n,1);      % hull polygon [x z]
TRI     = cell(n,1);      % full triangulation, reused by the depth sweep
LMK     = cell(n,1);      % geometry-derived landmarks

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

    tri  = stlread(char(stlPath(i)));
    TRI{i} = tri;
    LMK{i} = local_landmarks(tri, REF_LANDMARKS);
    tris = local_clip_below(tri, opt.CutY);
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

% ── depth sweep ──────────────────────────────────────────────────────────
SW = struct('z', [], 'ycut', [], 'Abare', [], 'Ahull', [], 'Alocal', []);
if opt.DepthSweep
    z  = opt.SweepZ(:);
    nz = numel(z);
    SW.z = z;  SW.ycut = LMK{1}.toeTip + z;
    SW.Abare = nan(nz,n); SW.Ahull = nan(nz,n); SW.Alocal = nan(nz,n);
    fprintf('\n--- depth sweep (%d cuts x %d models; polyshape unions, be patient) ---\n', nz, n);
    for i = 1:n
        t0 = tic;
        for k = 1:nz
            tris = local_clip_below(TRI{i}, SW.ycut(k));
            % At the shallowest cuts the plane only grazes the toe tips: the
            % clip is non-empty but every triangle is edge-on, so nothing is in
            % the bed yet and the answer is genuinely zero. local_footprint
            % must keep rejecting that case (there it means the cut plane or
            % the axis convention is wrong), so the sweep settles it here.
            if isempty(tris) || ~local_any_projected_area(tris)
                SW.Abare(k,i)=0; SW.Ahull(k,i)=0; SW.Alocal(k,i)=0; continue;
            end
            [~, SW.Abare(k,i)]  = local_footprint(tris);
            [~, SW.Ahull(k,i)]  = local_hull(tris);
            SW.Alocal(k,i)      = local_section(tris);
        end
        fprintf('  %-8s done in %.0f s\n', name(i), toc(t0));
    end
    % Validate against the independent sweep at the reference depths.
    fprintf('  validation vs trimesh/shapely 2026-08-22 (2%%):\n');
    okS = true;
    for j = 1:numel(REF_SWEEP_Z)
        kz = find(abs(z - REF_SWEEP_Z(j)) < 1e-6, 1);
        if isempty(kz), fprintf('    z = %g mm not on grid, skipped\n', REF_SWEEP_Z(j)); continue; end
        for i = 1:n
            [o1, m1] = local_check(sprintf('%s@%gmm', name(i), REF_SWEEP_Z(j)), 'A_bare',  SW.Abare(kz,i),  REF_SWEEP_BARE(j),   TOL_REL);
            [o2, m2] = local_check(sprintf('%s@%gmm', name(i), REF_SWEEP_Z(j)), 'A_hull',  SW.Ahull(kz,i),  REF_SWEEP_HULL(i,j), TOL_REL);
            [o3, m3] = local_check(sprintf('%s@%gmm', name(i), REF_SWEEP_Z(j)), 'a_local', SW.Alocal(kz,i), REF_SWEEP_LOC(j),    TOL_REL);
            fprintf('%s%s%s', m1, m2, m3);
            okS = okS && o1 && o2 && o3;
        end
    end
    if ~okS
        error('fig_foot_schematic:sweepValidationFailed', ...
              'Depth-sweep areas disagree with the reference by more than %g%%.', 100*TOL_REL);
    end
end

% ── figures ──────────────────────────────────────────────────────────────
if opt.Save && ~isfolder(opt.OutDir), mkdir(opt.OutDir); end

iDefault = find(name == "Default", 1);
figA = local_figure_assembly(char(stlPath(iDefault)), LMK{iDefault}, opt);
figB = local_figure_footprints(name, FP, HULL, Abare, Ahull, opt);
figC = gobjects(1);
if opt.DepthSweep
    figC = local_figure_sweep(name, SW, LMK{iDefault}, opt);
end

written = strings(0,1);
if opt.Save
    written = [written; local_export(figA, opt.OutDir, 'fig_foot_assembly')];
    written = [written; local_export(figB, opt.OutDir, 'fig_foot_footprints')];
    csvPath = local_write_csv(opt, name, stlPath, Abare, Ahull);
    written = [written; string(csvPath)];
    if opt.DepthSweep
        written = [written; local_export(figC, opt.OutDir, 'fig_foot_area_vs_depth')];
        written = [written; string(local_write_sweep_csv(opt, name, SW))];
    end
    fprintf('\n--- written ---\n');
    fprintf('  %s\n', written);
end

if ~opt.Show
    close(figA); close(figB);
    if isgraphics(figC), close(figC); end
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
R.landmarks      = LMK;
R.sweep          = SW;
R.written        = written;
R.figAssembly    = figA;
R.figFootprints  = figB;
R.figSweep       = figC;
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

function tf = local_any_projected_area(tris)
%LOCAL_ANY_PROJECTED_AREA  Does any clipped triangle project to area on XZ?
%   Separates "the cut only grazes the mesh" -- a real zero, at the toe tip --
%   from "the projection is degenerate", which local_footprint must keep
%   treating as a wrong cut plane. Same 1e-9 threshold local_footprint uses.
    x = tris(:, [1 4 7]);
    z = tris(:, [3 6 9]);
    a = 0.5 * abs( (x(:,2)-x(:,1)).*(z(:,3)-z(:,1)) ...
                 - (x(:,3)-x(:,1)).*(z(:,2)-z(:,1)) );
    tf = any(a >= 1e-9);
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

function areaCm2 = local_section(tris)
%LOCAL_SECTION  True cross-section area of the model at the cut plane.
%
%   Divergence theorem on the clipped OPEN surface: for the closed body below
%   the cut, the flux of the constant field (0,1,0) vanishes, so the flux
%   through the open surface equals minus the flux through the cap, and the
%   cap is exactly the section. Each triangle's flux is the y-component of its
%   vector area, 0.5*cross(B-A, C-A), which is its SIGNED projection onto XZ.
%   Requires a watertight, consistently outward-wound STL (true for the
%   OpenSCAD exports; the sign check below catches the opposite winding).
    A = tris(:,1:3); B = tris(:,4:6); C = tris(:,7:9);
    nv = 0.5 * cross(B - A, C - A, 2);
    areaCm2 = -sum(nv(:,2)) / 100;
    if areaCm2 < -1e-6
        error('fig_foot_schematic:winding', ...
              'Negative section area: STL faces are wound inward. Flip the sign in local_section or re-export.');
    end
    areaCm2 = max(areaCm2, 0);
end

function L = local_landmarks(tri, ref)
%LOCAL_LANDMARKS  Depth landmarks derived from the STL, asserted against the
%   documented values (cad/README.md). The selection rules depend on the STL
%   frame (drop = -Y, toes toward +x, mount at |z| = 10) and nothing else.
    P = tri.Points; x = P(:,1); y = P(:,2); z = P(:,3);
    L.toeTip      = min(y);
    L.barBottom   = min(y(x > 59 & x < 64));            % base block underside
    L.footTop     = max(y(x > 63 & y < -49));           % highest point of the toe structure
    L.postTop     = max(y(x > 58.5 & y < -10 & y > -49));
    L.mountBottom = min(y(abs(z) > 8));
    f = fieldnames(ref);
    for i = 1:numel(f)
        if abs(L.(f{i}) - ref.(f{i})) > 0.05
            error('fig_foot_schematic:landmark', ...
                  'Landmark %s = %.3f mm, documented %.3f mm. CAD revised? Update cad/README.md.', ...
                  f{i}, L.(f{i}), ref.(f{i}));
        end
    end
    L.z = @(yy) yy - L.toeTip;                           % toe-tip depth for a given y
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
function fig = local_figure_assembly(stlPath, L, opt)
%LOCAL_FIGURE_ASSEMBLY  Default model, FULL mesh, projected along Z onto XY.
%   Silhouette only: the union boundary as a black line over a light grey fill,
%   with no wireframe. Mesh lines would read as structure that is not there.

    ws = warning;                     % full current state, restored on exit
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    warning('off', 'MATLAB:polyshape:boundary3Points');
    restore = onCleanup(@() warning(ws)); %#ok<NASGU>

    tri = stlread(stlPath);
    P = tri.Points;  F = tri.ConnectivityList;

    parts = cell(size(F,1),1);
    kept = 0;
    for f = 1:size(F,1)
        x = P(F(f,:), 1).';
        y = P(F(f,:), 2).';
        if polyarea(x, y) < 1e-9, continue; end
        kept = kept + 1;
        parts{kept} = polyshape(x, y, 'Simplify', true);
    end
    pg = union([parts{1:kept}]);

    fig = figure('Color','w','Units','inches','Position',[1 1 3.4 5.0], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    ax = axes(fig); hold(ax,'on');

    plot(pg, 'Parent', ax, 'FaceColor', [0.85 0.85 0.85], 'FaceAlpha', 1, ...
             'EdgeColor', 'k', 'LineWidth', 1.0);

    % Equal aspect FIRST: it changes the limits, and every overlay below is
    % placed in data coordinates relative to them.
    axis(ax, 'equal');
    xl = xlim(ax);  yl = ylim(ax);
    xr = xl(2) - xl(1);

    % Cut plane, with its toe-tip depth so the reader can place it on z.
    plot(ax, xl, [opt.CutY opt.CutY], 'k--', 'LineWidth', 1.0);
    text(ax, xl(2), opt.CutY, sprintf(' area cut: foot/bar junction, z = %.2f cm', L.z(opt.CutY)/10), ...
         'FontSize', 8, 'VerticalAlignment', 'bottom', ...
         'HorizontalAlignment', 'right');
    % Geometry-derived depth landmarks (dotted), not typed-in positions.
    lmk = {'toe bars: bottom', L.barBottom; 'post: top', L.postTop};
    for i = 1:size(lmk,1)
        plot(ax, xl, [lmk{i,2} lmk{i,2}], ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
        text(ax, xl(1), lmk{i,2}, sprintf(' %s, z = %.2f cm', lmk{i,1}, L.z(lmk{i,2})/10), ...
             'FontSize', 7, 'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'bottom');
    end

    % Drop direction: an arrow along -Y, drawn in data coords so it cannot
    % drift from the geometry the way a figure-space annotation would.
    ax0 = xl(1) - 0.10*xr;
    yTop = yl(2) - 0.05*(yl(2)-yl(1));
    yBot = yl(1) + 0.20*(yl(2)-yl(1));
    quiver(ax, ax0, yTop, 0, yBot - yTop, 'Color', 'k', 'LineWidth', 1.2, ...
           'MaxHeadSize', 0.4, 'AutoScale', 'off');
    text(ax, ax0, (yTop+yBot)/2, ' drop direction', 'Rotation', 90, ...
         'FontSize', 8, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom');

    % 10 mm scale bar
    sbX = xl(2) - 0.05*xr - 10;
    sbY = yl(1) + 0.06*(yl(2)-yl(1));
    plot(ax, [sbX sbX+10], [sbY sbY], 'k-', 'LineWidth', 2.0);
    text(ax, sbX+5, sbY, '10 mm', 'FontSize', 8, ...
         'HorizontalAlignment','center', 'VerticalAlignment','top');

    if opt.Labels
        % Positions come from the geometry-derived landmarks, so the labels
        % name the part that is actually at that y. (The earlier fractional
        % placement put 'beam' beside the inclined bar, which is not the
        % structure the cut is defined on.)
        xr2 = xl(2) + 0.02*xr;
        lab = { 'mount',                         (L.mountBottom + yl(2))/2 ; ...
                'inclined bar (linkage)',        (opt.CutY + L.mountBottom)/2 ; ...
                'marker post',                   L.postTop ; ...
                'rectangular beam (proximal toe bars)', (opt.CutY + L.barBottom)/2 ; ...
                'distal toe segments',           (L.barBottom + L.toeTip)/2 };
        for i = 1:size(lab,1)
            text(ax, xr2, lab{i,2}, lab{i,1}, 'FontSize', 8, ...
                 'HorizontalAlignment','left', 'VerticalAlignment','middle');
        end
    end

    axis(ax, 'off');
    title(ax, 'Foot assembly (Default), side view', 'FontSize', 9);
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
    pad = 0.05 * max(max(xs)-min(xs), max(zs)-min(zs));
    XL = [min(xs)-pad, max(xs)+pad];
    ZL = [min(zs)-pad, max(zs)+pad];

    for i = 1:n
        ax = nexttile(tl); hold(ax,'on'); box(ax,'on');

        % Bare footprint. polyshape renders holes as background, so the gaps
        % between the toes read as gaps rather than as filled area.
        plot(FP{i}, 'Parent', ax, 'FaceColor', [0.75 0.75 0.75], ...
             'FaceAlpha', 1, 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);

        % Convex hull
        plot(ax, HULL{i}(:,1), HULL{i}(:,2), 'r--', 'LineWidth', 1.2);

        axis(ax, 'equal');
        xlim(ax, XL); ylim(ax, ZL);
        set(ax, 'FontSize', 8, 'LineWidth', 0.5, 'Layer', 'top');
        xlabel(ax, 'x (mm)', 'FontSize', 8);
        if i == 1, ylabel(ax, 'z (mm)', 'FontSize', 8); end

        % Title carries the COMPUTED values, never typed-in ones.
        title(ax, sprintf('%s: A_{bare} = %.2f, A_{hull} = %.2f cm^2', ...
                          name(i), Abare(i), Ahull(i)), ...
              'FontSize', 8, 'FontWeight','normal');
    end

    title(tl, sprintf('Intruding footprint at y = %g mm (grey) and convex hull (red)', ...
                      opt.CutY), 'FontSize', 9);
end

% ═════════════════════════════════════════════════════════════════════════
%  FIGURE C -- area in the bed vs toe-tip depth
% ═════════════════════════════════════════════════════════════════════════
function fig = local_figure_sweep(name, SW, L, opt)
%LOCAL_FIGURE_SWEEP  2x1: cumulative projected area (bare solid, hull dashed)
%   and local cross-section, vs toe-tip depth in cm. Landmarks dotted, cut
%   dashed. Everything below the surface is counted, so beyond the cut the
%   curves include the inclined bar and the post.
    n  = numel(name);
    zc = SW.z / 10;
    fig = figure('Color','w','Units','inches','Position',[1 1 3.4 4.6], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl = tiledlayout(fig, 2, 1, 'Padding','compact', 'TileSpacing','compact');
    col = lines(n);
    ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on');
    for i = 1:n
        plot(ax1, zc, SW.Ahull(:,i), '--', 'Color', col(i,:), 'LineWidth', 1.0, ...
             'DisplayName', sprintf('%s hull', name(i)));
        plot(ax1, zc, SW.Abare(:,i), '-',  'Color', col(i,:), 'LineWidth', 1.0, ...
             'DisplayName', sprintf('%s bare', name(i)));
    end
    ylabel(ax1, 'A below surface (cm^2)', 'FontSize', 8);
    legend(ax1, 'Location', 'northwest', 'FontSize', 6, 'NumColumns', 2);
    ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on');
    for i = 1:n
        plot(ax2, zc, SW.Alocal(:,i), '-', 'Color', col(i,:), 'LineWidth', 1.0);
    end
    ylabel(ax2, 'a(z) at surface (cm^2)', 'FontSize', 8);
    xlabel(ax2, 'toe-tip depth z (cm)', 'FontSize', 8);
    for ax = [ax1 ax2]
        yl = ylim(ax);
        for yy = [L.barBottom, L.postTop]
            plot(ax, L.z(yy)/10*[1 1], yl, ':', 'Color', [0.4 0.4 0.4]);
        end
        plot(ax, L.z(opt.CutY)/10*[1 1], yl, 'k--');
        ylim(ax, yl); xlim(ax, [0 max(zc)]);
        set(ax, 'FontSize', 8, 'LineWidth', 0.5, 'Layer', 'top');
    end
    title(tl, sprintf('Area in the bed vs depth (cut at z = %.2f cm dashed)', L.z(opt.CutY)/10), 'FontSize', 9);
end

function csvPath = local_write_sweep_csv(opt, name, SW)
    csvPath = fullfile(opt.OutDir, 'foot_area_vs_depth.csv');
    fid = fopen(csvPath, 'w');
    if fid < 0, error('fig_foot_schematic:csvFailed', 'Could not write %s', csvPath); end
    fprintf(fid, '# computed by scripts/fig_foot_schematic.m on %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# z_mm = toe-tip depth below the bed surface; y_cut_mm = STL cut plane\n');
    fprintf(fid, '# A_bare/A_hull: cumulative projected area of everything below the surface, cm^2\n');
    fprintf(fid, '# a_local: true cross-section at the surface plane, cm^2\n');
    fprintf(fid, 'z_mm,y_cut_mm');
    for i = 1:numel(name)
        fprintf(fid, ',A_bare_%s,A_hull_%s,a_local_%s', name(i), name(i), name(i));
    end
    fprintf(fid, '\n');
    for k = 1:numel(SW.z)
        fprintf(fid, '%.3f,%.3f', SW.z(k), SW.ycut(k));
        for i = 1:numel(name)
            fprintf(fid, ',%.5f,%.5f,%.5f', SW.Abare(k,i), SW.Ahull(k,i), SW.Alocal(k,i));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
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
