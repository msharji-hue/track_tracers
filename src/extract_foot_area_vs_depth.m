function out = extract_foot_area_vs_depth(stlFile, varargin)
% EXTRACT_FOOT_AREA_VS_DEPTH
%   Computes cross-sectional area A(z) for a foot penetrating in the -Y direction.
%   Slices the mesh at constant-y planes; cross-section lives in the x-z plane.
%
%   Usage:
%       out = extract_foot_area_vs_depth('foot.stl')
%       out = extract_foot_area_vs_depth('foot.stl', 'nLevels', 500, 'alpha', 1.5)
%       out = extract_foot_area_vs_depth('foot.stl', 'alpha', 1.5, 'gamma', 0.3)
%
%   Optional name-value args:
%       'nLevels'      - number of depth slices (default 300)
%       'alpha'        - fixed alphaShape radius in mm (default: auto)
%                        only used as fallback/visual reference now
%       'gamma'        - interstitial/jamming fraction for A_eff (default 0.3)
%       'sgWindow'     - Savitzky-Golay smoothing window (default 21, must be odd)
%       'sgOrder'      - SG polynomial order (default 3)
%       'verifyDepths' - vector of depths (cm) to plot cross-sections for
%                        visual verification (default [0.2 0.5 1.0 1.5])
%       'saveOutput'   - true/false to save out struct to .mat (default true)
%
%   Output fields:
%       depth_mm / depth_cm      - penetration depth from leading tip
%       A_bare_mm2 / cm2         - closed-contour material area
%       A_hull_mm2 / cm2         - convex hull envelope area
%       A_gap_mm2  / cm2         - interstitial gap area = hull - bare
%       A_eff_cm2                - effective area = bare + gamma*gap
%       A_bare_sm / A_hull_sm    - smoothed area curves
%       A_gap_sm / A_eff_sm      - smoothed gap/effective area curves

    %% ── Parse inputs ──────────────────────────────────────────────────────
    p = inputParser;
    addParameter(p, 'nLevels',      300,               @isnumeric);
    addParameter(p, 'alpha',        [],                @(x) isempty(x)||isnumeric(x));
    addParameter(p, 'gamma',        0.3,               @isnumeric);
    addParameter(p, 'sgWindow',     21,                @isnumeric);
    addParameter(p, 'sgOrder',      3,                 @isnumeric);
    addParameter(p, 'verifyDepths', [0.2 0.5 1.0 1.5], @isnumeric);
    addParameter(p, 'saveOutput',   true,              @islogical);
    parse(p, varargin{:});
    opt = p.Results;

    gamma = opt.gamma;

    %% ── Load STL ──────────────────────────────────────────────────────────
    fprintf('Loading STL: %s\n', stlFile);
    TR = stlread(stlFile);
    F  = TR.ConnectivityList;
    V  = TR.Points;
    fprintf('  %d faces,  %d vertices\n', size(F,1), size(V,1));

    % Penetration axis is -Y; cross-section lives in x-z plane
    y_tip = min(V(:,2));   % leading edge / first contact
    y_top = max(V(:,2));   % trailing edge

    fprintf('  Y range: %.2f to %.2f mm  (%.2f mm total)\n', ...
        y_tip, y_top, y_top - y_tip);

    %% ── Depth levels ──────────────────────────────────────────────────────
    nL       = opt.nLevels;
    yLevels  = linspace(y_tip, y_top, nL);
    depth_mm = (yLevels - y_tip)';
    depth_cm = depth_mm / 10;

    %% ── Precompute face y-ranges for fast filtering ───────────────────────
    y_verts = V(:,2);
    y_face  = y_verts(F);
    y_fmin  = min(y_face, [], 2);
    y_fmax  = max(y_face, [], 2);

    % Precompute vertex coordinates for all face vertices
    v1 = V(F(:,1), :);
    v2 = V(F(:,2), :);
    v3 = V(F(:,3), :);

    %% ── Auto alpha if not specified ───────────────────────────────────────
    foot_width = max(V(:,1)) - min(V(:,1));
    alpha_auto = foot_width * 0.08;

    if isempty(opt.alpha)
        alpha_use = alpha_auto;
        fprintf('  Auto alpha = %.2f mm  (fallback only)\n', alpha_use);
    else
        alpha_use = opt.alpha;
        fprintf('  Using fixed alpha = %.2f mm  (fallback only)\n', alpha_use);
    end

    fprintf('  Bare area method: closed slice contours\n');
    fprintf('  Using gamma = %.2f for A_eff = A_bare + gamma*A_gap\n', gamma);

    %% ── Main slicing loop ─────────────────────────────────────────────────
    A_bare_mm2 = nan(nL, 1);
    A_hull_mm2 = nan(nL, 1);
    nLoops     = zeros(nL, 1);

    fprintf('Slicing %d depth levels...\n', nL);
    t_start = tic;

    for k = 1:nL
        y0 = yLevels(k);

        % Only use faces whose y-range straddles the slicing plane
        straddle = y_fmin <= y0 & y_fmax >= y0;
        if ~any(straddle)
            continue;
        end

        sv1 = v1(straddle,:);
        sv2 = v2(straddle,:);
        sv3 = v3(straddle,:);

        % Get actual line segments from triangle-plane intersections
        [segA, segB, pts] = triangle_slice_segments(sv1, sv2, sv3, y0);

        if size(pts,1) < 3
            continue;
        end

        pts = unique(round(pts, 5), 'rows');

        if size(pts,1) < 3
            continue;
        end

        % Convex hull area: occupied envelope
        try
            K = convhull(pts(:,1), pts(:,2));
            A_hull_mm2(k) = polyarea(pts(K,1), pts(K,2));
        catch
            A_hull_mm2(k) = NaN;
        end

        % Bare material area: sum of closed cross-section loops
        try
            [A_bare_mm2(k), loops] = area_from_closed_segments(segA, segB);
            nLoops(k) = numel(loops);
        catch
            A_bare_mm2(k) = NaN;
        end

        % Fallback only if closed-loop reconstruction fails
        if ~isfinite(A_bare_mm2(k)) || A_bare_mm2(k) <= 0
            try
                shp = alphaShape(pts(:,1), pts(:,2), alpha_use);
                A_bare_mm2(k) = area(shp);
            catch
                A_bare_mm2(k) = NaN;
            end
        end

        % Physical cleanup at the raw-area level
        if isfinite(A_bare_mm2(k))
            A_bare_mm2(k) = max(A_bare_mm2(k), 0);
        end

        if isfinite(A_hull_mm2(k))
            A_hull_mm2(k) = max(A_hull_mm2(k), 0);
        end

        if isfinite(A_bare_mm2(k)) && isfinite(A_hull_mm2(k))
            A_bare_mm2(k) = min(A_bare_mm2(k), A_hull_mm2(k));
        end

        if mod(k,50)==0 || k==nL
            fprintf('  %d / %d  (%.1f s elapsed)\n', k, nL, toc(t_start));
        end
    end

    A_gap_mm2 = max(A_hull_mm2 - A_bare_mm2, 0);
    A_eff_mm2 = A_bare_mm2 + gamma*A_gap_mm2;

    %% ── Convert to cm² ───────────────────────────────────────────────────
    A_bare_cm2 = A_bare_mm2 / 100;
    A_hull_cm2 = A_hull_mm2 / 100;
    A_gap_cm2  = A_gap_mm2  / 100;
    A_eff_cm2  = A_eff_mm2  / 100;

    %% ── Smooth with Savitzky-Golay filter ────────────────────────────────
    sgW = opt.sgWindow;
    sgO = opt.sgOrder;

    if mod(sgW,2)==0
        sgW = sgW + 1;
    end

    A_bare_sm = nan(nL,1);
    A_hull_sm = nan(nL,1);

    ok_bare = isfinite(A_bare_cm2);
    ok_hull = isfinite(A_hull_cm2);

    if sum(ok_bare) > sgW
        A_bare_sm(ok_bare) = sgolayfilt(A_bare_cm2(ok_bare), sgO, sgW);
    end

    if sum(ok_hull) > sgW
        A_hull_sm(ok_hull) = sgolayfilt(A_hull_cm2(ok_hull), sgO, sgW);
    end

    % Physical cleanup after smoothing
    A_bare_sm = max(A_bare_sm, 0);
    A_hull_sm = max(A_hull_sm, 0);

    bad = A_bare_sm > A_hull_sm;
    A_bare_sm(bad) = A_hull_sm(bad);

    A_gap_sm = max(A_hull_sm - A_bare_sm, 0);
    A_eff_sm = A_bare_sm + gamma*A_gap_sm;

    %% ── Package output ───────────────────────────────────────────────────
    out.depth_mm   = depth_mm;
    out.depth_cm   = depth_cm;

    out.A_bare_mm2 = A_bare_mm2;
    out.A_hull_mm2 = A_hull_mm2;
    out.A_gap_mm2  = A_gap_mm2;
    out.A_eff_mm2  = A_eff_mm2;

    out.A_bare_cm2 = A_bare_cm2;
    out.A_hull_cm2 = A_hull_cm2;
    out.A_gap_cm2  = A_gap_cm2;
    out.A_eff_cm2  = A_eff_cm2;

    out.A_bare_sm  = A_bare_sm;
    out.A_hull_sm  = A_hull_sm;
    out.A_gap_sm   = A_gap_sm;
    out.A_eff_sm   = A_eff_sm;

    out.nLoops     = nLoops;
    out.alpha_used = alpha_use;
    out.gamma      = gamma;
    out.y_tip      = y_tip;
    out.y_top      = y_top;
    out.bare_area_method = 'closed slice contours';

    %% ── Console summary ──────────────────────────────────────────────────
    [peakBare, idxBare] = max(A_bare_sm, [], 'omitnan');
    [peakHull, idxHull] = max(A_hull_sm, [], 'omitnan');
    [peakGap,  idxGap]  = max(A_gap_sm,  [], 'omitnan');
    [peakEff,  idxEff]  = max(A_eff_sm,  [], 'omitnan');

    fprintf('\nDone. Total time: %.1f s\n', toc(t_start));
    fprintf('Depth range:      0 to %.2f cm\n', max(depth_cm));
    fprintf('Peak bare area:   %.4f cm^2  at z = %.2f cm\n', ...
        peakBare, depth_cm(idxBare));
    fprintf('Peak hull area:   %.4f cm^2  at z = %.2f cm\n', ...
        peakHull, depth_cm(idxHull));
    fprintf('Peak gap area:    %.4f cm^2  at z = %.2f cm\n', ...
        peakGap, depth_cm(idxGap));
    fprintf('Peak eff. area:   %.4f cm^2  at z = %.2f cm  (gamma = %.2f)\n', ...
        peakEff, depth_cm(idxEff), gamma);

    %% ── Main area plot ────────────────────────────────────────────────────
    fig1 = figure('Name','A(z) vs penetration depth', ...
                  'ToolBar','none','MenuBar','none');
    fig1.Position = [100 100 700 480];

    ax = axes(fig1, 'Position', [0.13 0.13 0.78 0.80]);
    hold(ax,'on');

    % Raw faint curves
    plot(ax, depth_cm, A_bare_cm2, '-', 'Color', [0.3 0.6 0.9 0.3], ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(ax, depth_cm, A_hull_cm2, '-', 'Color', [0.9 0.4 0.1 0.3], ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');

    % Smoothed bold curves
    plot(ax, depth_cm, A_bare_sm, '-',  'Color', [0.15 0.45 0.80], ...
        'LineWidth', 2.2, 'DisplayName', 'Bare foot area');
    plot(ax, depth_cm, A_hull_sm, '--', 'Color', [0.85 0.30 0.05], ...
        'LineWidth', 2.2, 'DisplayName', 'Convex hull area');
    plot(ax, depth_cm, A_gap_sm,  ':',  'Color', [0.20 0.65 0.30], ...
        'LineWidth', 2.0, 'DisplayName', 'Interstitial gap');
    plot(ax, depth_cm, A_eff_sm,  '-.', 'Color', [0.20 0.20 0.20], ...
        'LineWidth', 2.0, 'DisplayName', sprintf('Effective area, \\gamma = %.2f', gamma));

    set(ax, 'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in');

    grid(ax,'off');

    xlabel(ax, 'Penetration depth,  $z$  (cm)', ...
        'FontSize',15,'Interpreter','latex','Color',[0 0 0]);

    ylabel(ax, 'Cross-sectional area,  $A$  (cm$^2$)', ...
        'FontSize',15,'Interpreter','latex','Color',[0 0 0]);

    legend(ax,'show','FontSize',11,'Box','on','Location','best', ...
        'Interpreter','tex','EdgeColor',[0.25 0.25 0.25]);

    %% ── Verification: cross-sections at selected depths ──────────────────
    vd = opt.verifyDepths;
    nV = numel(vd);

    if nV > 0
        fig2 = figure('Name','Cross-section verification', ...
                      'ToolBar','none','MenuBar','none');
        fig2.Position = [150 150 280*nV 260];

        for vi = 1:nV
            target_cm = vd(vi);
            y0 = y_tip + target_cm*10;

            straddle = y_fmin <= y0 & y_fmax >= y0;

            sv1 = v1(straddle,:);
            sv2 = v2(straddle,:);
            sv3 = v3(straddle,:);

            [segA, segB, pts] = triangle_slice_segments(sv1, sv2, sv3, y0);

            axv = subplot(1, nV, vi);
            hold(axv,'on');

            if size(pts,1) >= 3
                pts = unique(round(pts,5),'rows');

                % Hull fill
                try
                    K = convhull(pts(:,1), pts(:,2));
                    fill(axv, pts(K,1), pts(K,2), [0.95 0.85 0.70], ...
                        'EdgeColor', [0.85 0.30 0.05], ...
                        'LineWidth', 1.2, ...
                        'FaceAlpha', 0.4);
                catch
                end

                % Closed-loop bare area fill
                try
                    [~, loops] = area_from_closed_segments(segA, segB);

                    for li = 1:numel(loops)
                        L = loops{li};
                        if size(L,1) >= 3
                            fill(axv, L(:,1), L(:,2), [0.15 0.45 0.80], ...
                                'FaceAlpha', 0.6, ...
                                'EdgeColor', [0.10 0.30 0.60], ...
                                'LineWidth', 1.0);
                        end
                    end
                catch
                    % Fallback visual only
                    try
                        shp = alphaShape(pts(:,1), pts(:,2), alpha_use);
                        plot(shp, 'Parent', axv, ...
                            'FaceColor', [0.15 0.45 0.80], ...
                            'FaceAlpha', 0.6, ...
                            'EdgeColor', [0.10 0.30 0.60], ...
                            'LineWidth', 1.0);
                    catch
                    end
                end

                plot(axv, pts(:,1), pts(:,2), 'k.', 'MarkerSize', 2);
            end

            axis(axv,'equal');

            set(axv,'FontSize',10,'Box','on','TickDir','in');

            title(axv, sprintf('z = %.1f cm', target_cm), ...
                'FontSize',10,'Interpreter','tex');

            xlabel(axv,'x (mm)','FontSize',9);

            if vi == 1
                ylabel(axv,'z_{\rm STL} (mm)','FontSize',9);
            end
        end

        sgtitle('Cross-section verification  (blue = closed-loop bare, tan = hull)', ...
            'FontSize', 11);
    end

    %% ── Save ─────────────────────────────────────────────────────────────
    if opt.saveOutput
        [fPath, fName] = fileparts(stlFile);

        if isempty(fPath)
            fPath = pwd;
        end

        outFile = fullfile(fPath, [fName '_area_profile.mat']);
        save(outFile, 'out', '-v7.3');
        fprintf('Saved: %s\n', outFile);

        pdfFile = fullfile(fPath, [fName '_area_profile.pdf']);
        exportgraphics(fig1, pdfFile, ...
            'ContentType','vector','BackgroundColor','white');
        fprintf('Saved: %s\n', pdfFile);
    end
end

%% ── Helper: triangle-plane intersection segments ─────────────────────────
function [segA, segB, ptsAll] = triangle_slice_segments(va, vb, vc, y0)
% For each triangle, find the segment made by intersecting the triangle
% with plane y = y0.
%
% Returned coordinates are [x, z_STL].

    nT = size(va,1);

    segA = [];
    segB = [];
    ptsAll = [];

    for ii = 1:nT
        tri = [va(ii,:); vb(ii,:); vc(ii,:)];

        pLocal = [];

        edgePairs = [1 2; 2 3; 3 1];

        for ee = 1:3
            i1 = edgePairs(ee,1);
            i2 = edgePairs(ee,2);

            p1 = tri(i1,:);
            p2 = tri(i2,:);

            y1 = p1(2);
            y2 = p2(2);

            if y1 == y2
                continue;
            end

            crosses = (y1 - y0) * (y2 - y0) <= 0;

            if crosses
                t = (y0 - y1) / (y2 - y1);

                if t >= -1e-10 && t <= 1 + 1e-10
                    p = p1 + t*(p2 - p1);
                    pLocal = [pLocal; p(1), p(3)];
                end
            end
        end

        if size(pLocal,1) >= 2
            pLocal = unique(round(pLocal, 6), 'rows');

            if size(pLocal,1) >= 2
                % If more than two points occur because the plane passes
                % through a vertex, keep the farthest pair as the segment.
                D = squareform(pdist(pLocal));
                [~, imax] = max(D(:));
                [r, c] = ind2sub(size(D), imax);

                pa = pLocal(r,:);
                pb = pLocal(c,:);

                if norm(pa - pb) > 1e-8
                    segA = [segA; pa];
                    segB = [segB; pb];
                    ptsAll = [ptsAll; pa; pb];
                end
            end
        end
    end
end

%% ── Helper: closed-loop area from unordered slice segments ───────────────
function [A_total, loops] = area_from_closed_segments(segA, segB)
% Reconstruct closed loops from unordered slice segments and sum loop areas.
% This is intended for cross-sections of a watertight STL.
%
% segA, segB are [nSegments x 2] endpoints in the x-z plane.

    A_total = NaN;
    loops = {};

    if isempty(segA) || isempty(segB)
        return;
    end

    tol = 1e-5;

    P = [segA; segB];
    Pr = round(P / tol) * tol;

    [nodes, ~, ic] = unique(Pr, 'rows', 'stable');

    nSeg = size(segA,1);
    e1 = ic(1:nSeg);
    e2 = ic(nSeg+1:end);

    edges = [e1, e2];

    % Remove zero-length edges
    edges(edges(:,1) == edges(:,2), :) = [];

    if isempty(edges)
        return;
    end

    % Remove duplicate undirected edges
    edgesSorted = sort(edges, 2);
    [~, ia] = unique(edgesSorted, 'rows', 'stable');
    edges = edges(ia,:);

    nEdges = size(edges,1);
    used = false(nEdges,1);

    for startEdge = 1:nEdges

        if used(startEdge)
            continue;
        end

        a = edges(startEdge,1);
        b = edges(startEdge,2);

        loopNodes = [a; b];
        used(startEdge) = true;

        current = b;
        previous = a;

        closed = false;

        while true
            % Candidate unused edges attached to current node
            cand = find(~used & (edges(:,1) == current | edges(:,2) == current));

            if isempty(cand)
                break;
            end

            % Prefer edge that does not immediately return to previous,
            % unless that is the only available option.
            nextEdge = cand(1);
            nextNode = other_node(edges(nextEdge,:), current);

            if numel(cand) > 1
                for cc = 1:numel(cand)
                    testNode = other_node(edges(cand(cc),:), current);
                    if testNode ~= previous
                        nextEdge = cand(cc);
                        nextNode = testNode;
                        break;
                    end
                end
            end

            used(nextEdge) = true;

            if nextNode == loopNodes(1)
                closed = true;
                break;
            end

            loopNodes = [loopNodes; nextNode];

            previous = current;
            current = nextNode;

            if numel(loopNodes) > nEdges + 2
                break;
            end
        end

        if closed && numel(loopNodes) >= 3
            xy = nodes(loopNodes, :);

            % Remove repeated consecutive points if any
            dxy = [true; any(abs(diff(xy,1,1)) > tol, 2)];
            xy = xy(dxy,:);

            if size(xy,1) >= 3
                loops{end+1} = xy; %#ok<AGROW>
            end
        end
    end

    if isempty(loops)
        A_total = NaN;
        return;
    end

    A_total = 0;

    for ii = 1:numel(loops)
        xy = loops{ii};

        if size(xy,1) >= 3
            A_total = A_total + abs(polyarea(xy(:,1), xy(:,2)));
        end
    end
end

%% ── Helper: get opposite node on an edge ─────────────────────────────────
function b = other_node(edge, a)
    if edge(1) == a
        b = edge(2);
    else
        b = edge(1);
    end
end