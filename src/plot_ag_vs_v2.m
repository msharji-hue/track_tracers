function [fig, d1_each, Fz_over_m, fitStats] = plot_ag_vs_v2(heights, ~, z_targets)
% PLOT_AG_VS_V2  Case 1b: raw Katsuragi-style diagnostic.
%
%   Main plot : a+g vs v
%   Inset     : a+g vs v^2
%
%   Case 1b:
%       a + g = F(z)/m + v^2/d1(z)
%
%   Each depth gets:
%       - its own intercept F(z)/m
%       - its own inertial slope 1/d1(z)
%
%   This is a diagnostic fit. It tests whether each fixed-depth curve is
%   v^2-like and whether the fitted d1 values are similar across depths.
%
%   Outputs:
%       fig        - figure handle
%       d1_each    - depth-specific inertial length scales [cm]
%       Fz_over_m  - fitted depth-dependent force terms [cm s^-2]
%       fitStats   - RMSE, R^2, residuals, prediction, etc.
%
%   Usage:
%       z_targets = [1.30 1.70 2.00];
%       [fig, d1_each, Fz, stats] = plot_ag_vs_v2(heights, cmap, z_targets);

    %% ── Settings ─────────────────────────────────────────────────────────
    nZ     = numel(z_targets);
    tol    = 0.05;     % depth tolerance [cm]
    ag_cap = 12000;    % acceleration cap [cm/s^2]
    nBins  = 16;

    %% ── Depth-level colors ───────────────────────────────────────────────
    depth_cols = [0.05 0.05 0.05;
                  0.55 0.10 0.25;
                  0.88 0.25 0.55;
                  0.98 0.65 0.80];

    if nZ > 4
        depth_cols = interp1(linspace(0,1,4), depth_cols, linspace(0,1,nZ));
    end

    %% ── Collect v and a+g at each fixed depth ────────────────────────────
    v_at_z  = cell(nZ,1);
    ag_at_z = cell(nZ,1);

    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials

            k   = heights(j).trials(i).kinematics;
            idx = k.impact_index:k.stopFrame;

            z  = k.z_smooth(idx);
            v  = k.v_smooth(idx);
            ag = k.a_plus_g(idx);

            for zi = 1:nZ

                in = abs(z - z_targets(zi)) < tol & ...
                     isfinite(v) & isfinite(ag) & ...
                     v > 0 & ag > 0 & ag < ag_cap;

                if any(in)
                    v_at_z{zi}  = [v_at_z{zi};  v(in)];
                    ag_at_z{zi} = [ag_at_z{zi}; ag(in)];
                end
            end
        end
    end

    %% ── Case 1b fit: one intercept and one v² slope per depth ────────────
    n_total = sum(cellfun(@numel, ag_at_z));

    if n_total == 0
        error('No valid data found at the requested z_targets.');
    end

    X_fit = zeros(n_total, 2*nZ);
    b_vec = zeros(n_total, 1);

    row = 1;

    for zi = 1:nZ
        n_zi = numel(ag_at_z{zi});

        if n_zi == 0
            continue;
        end

        rows = row:row+n_zi-1;

        X_fit(rows, zi)      = 1;              % depth-specific F(z)/m
        X_fit(rows, nZ + zi) = v_at_z{zi}.^2;  % depth-specific v^2 slope

        b_vec(rows) = ag_at_z{zi};

        row = row + n_zi;
    end

    params = lsqnonneg(X_fit, b_vec);

    Fz_over_m = params(1:nZ);
    slope_v2  = params(nZ+1:end);
    d1_each   = 1 ./ max(slope_v2, 1e-10);

    pred  = X_fit * params;
    resid = b_vec - pred;

    rmse   = sqrt(mean(resid.^2, 'omitnan'));
    ss_res = sum(resid.^2, 'omitnan');
    ss_tot = sum((b_vec - mean(b_vec,'omitnan')).^2, 'omitnan');
    r2     = 1 - ss_res/ss_tot;

    fitStats.rmse      = rmse;
    fitStats.r2        = r2;
    fitStats.resid     = resid;
    fitStats.pred      = pred;
    fitStats.b         = b_vec;
    fitStats.params    = params;
    fitStats.slope_v2  = slope_v2;
    fitStats.d1_each   = d1_each;
    fitStats.Fz_over_m = Fz_over_m;
    fitStats.z_targets = z_targets(:);

    fprintf('\n── Case 1b: Raw fit with depth-specific slopes ───────────────\n');
    fprintf('Model: a+g = F(z)/m + v^2/d1(z)\n');
    fprintf('RMSE = %.2f cm s^{-2}\n', rmse);
    fprintf('R^2  = %.4f\n', r2);

    for zi = 1:nZ
        fprintf('  z=%.2f cm: F/m=%.1f cm s^{-2}, d1=%.3f cm\n', ...
            z_targets(zi), Fz_over_m(zi), d1_each(zi));
    end

    fprintf('──────────────────────────────────────────────────────────────\n\n');

    %% ── Bin raw data ─────────────────────────────────────────────────────
    nonempty = ~cellfun(@isempty, v_at_z);

    if ~any(nonempty)
        error('No valid nonempty depth groups found.');
    end

    v_global_max = max(cellfun(@max, v_at_z(nonempty)));
    v_edges      = linspace(0, v_global_max*1.02, nBins+1);
    v_centers    = (v_edges(1:end-1) + v_edges(2:end)) / 2;

    bin_mean = nan(nZ, nBins);
    bin_std  = nan(nZ, nBins);

    for zi = 1:nZ
        if isempty(ag_at_z{zi})
            continue;
        end

        for b = 1:nBins
            in_b = v_at_z{zi} >= v_edges(b) & v_at_z{zi} < v_edges(b+1);

            if sum(in_b) >= 3
                bin_mean(zi,b) = mean(ag_at_z{zi}(in_b));
                bin_std(zi,b)  = std( ag_at_z{zi}(in_b));
            end
        end
    end

    ag_max_plot = max(bin_mean(:), [], 'omitnan') * 1.20;

    %% ── Main figure: a+g vs v ────────────────────────────────────────────
    fig = figure('Name','Case 1b: raw a+g vs v', ...
                 'ToolBar','none','MenuBar','none');

    fig.Position = [100 100 740 560];

    ax = axes(fig, 'Position', [0.11 0.12 0.82 0.83]);
    hold(ax,'on');

    v_fit = linspace(0, v_global_max*1.04, 400);

    for zi = 1:nZ
        if isempty(ag_at_z{zi})
            continue;
        end

        col = depth_cols(zi,:);

        % Faint raw scatter
        scatter(ax, v_at_z{zi}, ag_at_z{zi}, 5, col, 'filled', ...
            'MarkerFaceAlpha', 0.10, ...
            'MarkerEdgeColor', 'none', ...
            'HandleVisibility','off');

        % Binned means ± std
        ok = isfinite(bin_mean(zi,:));

        errorbar(ax, v_centers(ok), bin_mean(zi,ok), bin_std(zi,ok), ...
            'o', ...
            'Color', col, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'w', ...
            'MarkerSize', 9, ...
            'LineWidth', 1.4, ...
            'CapSize', 3, ...
            'DisplayName', sprintf('$z = %.2f$ cm', z_targets(zi)));

        % Depth-specific dashed fit line
        ag_fit = Fz_over_m(zi) + v_fit.^2 ./ d1_each(zi);

        plot(ax, v_fit, ag_fit, 'k--', ...
            'LineWidth', 1.8, ...
            'HandleVisibility','off');
    end

    set(ax,'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0], ...
        'YColor',[0 0 0], ...
        'XMinorTick','on', ...
        'YMinorTick','on', ...
        'TickDir','in', ...
        'XLim',[0, v_global_max*1.05], ...
        'YLim',[0, ag_max_plot]);

    grid(ax,'off');

    xlabel(ax,'$v$  (cm s$^{-1}$)', ...
        'FontSize',16, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);

    ylabel(ax,'$a+g$  (cm s$^{-2}$)', ...
        'FontSize',16, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);

    lgd = legend(ax,'show', ...
        'FontSize',11, ...
        'Box','on', ...
        'Interpreter','latex', ...
        'EdgeColor',[0.25 0.25 0.25], ...
        'Color',[1 1 1]);

    lgd.Location = 'southeast';

    d1_text = sprintf('%.2f, ', d1_each);
    d1_text = d1_text(1:end-2);

    text(ax, 0.97, 0.05, ...
        sprintf('RMSE = %.0f cm s$^{-2}$\n$R^2 = %.3f$\n$d_1(z)$ = [%s] cm', ...
        rmse, r2, d1_text), ...
        'FontSize',10, ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom', ...
        'FontAngle','italic', ...
        'Interpreter','latex', ...
        'Color',[0.15 0.15 0.15], ...
        'BackgroundColor',[1 1 1], ...
        'EdgeColor',[0.25 0.25 0.25], ...
        'LineWidth',1.0, ...
        'Margin',4);

    %% ── Inset: a+g vs v² ─────────────────────────────────────────────────
    ax_in = axes(fig, 'Position', [0.14 0.57 0.25 0.32]);
    hold(ax_in,'on');

    v2_max_all = max(cellfun(@(v) max(v.^2), v_at_z(nonempty)));
    v2_edges   = linspace(0, v2_max_all*1.02, nBins+1);
    v2_centers = (v2_edges(1:end-1) + v2_edges(2:end)) / 2;

    for zi = 1:nZ
        if isempty(ag_at_z{zi})
            continue;
        end

        col = depth_cols(zi,:);
        v2  = v_at_z{zi}.^2;

        bm2 = nan(1,nBins);

        for b = 1:nBins
            in_b = v2 >= v2_edges(b) & v2 < v2_edges(b+1);

            if sum(in_b) >= 3
                bm2(b) = mean(ag_at_z{zi}(in_b));
            end
        end

        ok2 = isfinite(bm2);

        plot(ax_in, v2_centers(ok2), bm2(ok2), 'o', ...
            'Color', col, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerSize', 6, ...
            'HandleVisibility','off');

        v2_fit = linspace(0, v2_max_all*1.04, 300);
        ag_fit = Fz_over_m(zi) + v2_fit ./ d1_each(zi);

        plot(ax_in, v2_fit, ag_fit, 'k--', ...
            'LineWidth', 1.4, ...
            'HandleVisibility','off');
    end

    set(ax_in,'FontSize',10,'Box','on','LineWidth',1.0, ...
        'XColor',[0 0 0], ...
        'YColor',[0 0 0], ...
        'TickDir','in', ...
        'XMinorTick','on', ...
        'YMinorTick','on', ...
        'XLim',[0, v2_max_all*1.04], ...
        'YLim',[0, max(bin_mean(:),[],'omitnan')*1.15]);

    grid(ax_in,'off');

    xlabel(ax_in,'$v^2$  (cm$^2$ s$^{-2}$)', ...
        'FontSize',10, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);

    ylabel(ax_in,'$a+g$  (cm s$^{-2}$)', ...
        'FontSize',10, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);
end