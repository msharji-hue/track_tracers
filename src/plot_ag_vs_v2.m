function [fig, d1_fit, Fz_over_m, fitStats] = plot_ag_vs_v2(heights, ~, z_targets, v_min)
% PLOT_AG_VS_V2  Katsuragi & Durian (2007) Fig. 3a — exact approach.
%
%   Fit: a+g = F(zi)/m + v²/d1
%       - SHARED d1 across ALL depths (same slope)
%       - Depth-specific intercepts F(zi)/m
%       - Parallel lines in v² inset = proof d1 is constant
%
%   v_min clips both fit and display. Excluded points not shown.
%   x-axis uses 50 cm/s tick spacing so the 0→v_min empty zone
%   is a minor fraction of the first interval rather than a full one.
%
%   Inputs:
%       heights    - struct array from group_trials_by_height
%       ~          - colormap (unused)
%       z_targets  - fixed depths to sample (cm), e.g. [1.0 1.4 1.7 2.1]
%       v_min      - minimum velocity for fit AND display (cm/s), default 0
%
%   Outputs:
%       fig        - figure handle
%       d1_fit     - shared inertial length scale [cm]
%       Fz_over_m  - fitted depth-specific intercepts [cm s^-2]
%       fitStats   - includes Fz_over_m_se for vertical error bars in Fig. 3b

    if nargin < 4
        v_min = 0;
    end

    nZ     = numel(z_targets);
    tol    = 0.05;
    ag_cap = 12000;
    nBins  = 16;

    %% ── Depth-level colors ───────────────────────────────────────────────
    depth_cols = [0.05 0.05 0.05;
                  0.55 0.10 0.25;
                  0.88 0.25 0.55;
                  0.98 0.65 0.80];

    if nZ > 4
        depth_cols = interp1(linspace(0,1,4), depth_cols, linspace(0,1,nZ));
    end

    %% ── Collect (v, a+g) at each fixed depth ─────────────────────────────
    v_raw  = cell(nZ,1);
    ag_raw = cell(nZ,1);

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
                    v_raw{zi}  = [v_raw{zi};  v(in)];
                    ag_raw{zi} = [ag_raw{zi}; ag(in)];
                end
            end
        end
    end

    %% ── Apply v_min to BOTH fit and display ──────────────────────────────
    v_at_z  = cell(nZ,1);
    ag_at_z = cell(nZ,1);

    for zi = 1:nZ
        ok = v_raw{zi} >= v_min;

        v_at_z{zi}  = v_raw{zi}(ok);
        ag_at_z{zi} = ag_raw{zi}(ok);
    end

    fprintf('v_min = %.1f cm/s  (applied to fit and display)\n', v_min);

    for zi = 1:nZ
        fprintf('  z=%.2f cm: %d total -> %d displayed/fitted\n', ...
            z_targets(zi), numel(v_raw{zi}), numel(v_at_z{zi}));
    end

    %% ── Katsuragi fit: SHARED d1, depth-specific F(zi)/m ─────────────────
    n_fit = sum(cellfun(@numel, ag_at_z));

    if n_fit == 0
        error('No valid data after v_min cutoff. Lower v_min.');
    end

    A_mat = zeros(n_fit, nZ+1);
    b_vec = zeros(n_fit, 1);
    row   = 1;

    for zi = 1:nZ

        n_zi = numel(ag_at_z{zi});

        if n_zi == 0
            continue;
        end

        rows = row:row+n_zi-1;

        A_mat(rows, zi)   = 1;              % depth-specific intercept
        A_mat(rows, nZ+1) = v_at_z{zi}.^2;  % shared v² slope

        b_vec(rows) = ag_at_z{zi};

        row = row + n_zi;
    end

    params    = lsqnonneg(A_mat, b_vec);
    Fz_over_m = params(1:nZ);
    d1_fit    = 1 / max(params(nZ+1), 1e-10);

    pred_fit  = A_mat * params;
    resid_fit = b_vec - pred_fit;

    rmse   = sqrt(mean(resid_fit.^2));
    ss_res = sum(resid_fit.^2);
    ss_tot = sum((b_vec - mean(b_vec)).^2);
    r2     = 1 - ss_res/ss_tot;

    %% ── Exact regression-based uncertainty for F(z_i)/m intercepts ──────
    % Model: a+g = F(z_i)/m + v²/d1
    % The first nZ parameters are the intercepts F(z_i)/m.
    n_params = size(A_mat, 2);
    dof      = n_fit - n_params;

    Fz_over_m_se = nan(nZ,1);

    if dof > 0
        sigma2     = sum(resid_fit.^2) / dof;
        cov_params = sigma2 * pinv(A_mat' * A_mat);
        param_se   = sqrt(diag(cov_params));

        Fz_over_m_se = param_se(1:nZ);
    end

    %% ── Package fit stats ────────────────────────────────────────────────
    fitStats.d1_fit       = d1_fit;
    fitStats.Fz_over_m    = Fz_over_m;
    fitStats.Fz_over_m_se = Fz_over_m_se;
    fitStats.rmse         = rmse;
    fitStats.r2           = r2;
    fitStats.v_min        = v_min;
    fitStats.z_targets    = z_targets(:);
    fitStats.params       = params;
    fitStats.pred_fit     = pred_fit;
    fitStats.resid_fit    = resid_fit;
    fitStats.A_mat        = A_mat;
    fitStats.b_vec        = b_vec;

    fprintf('\n-- Katsuragi fit: shared d1, depth-specific F(zi)/m ---------\n');
    fprintf('d1   = %.3f cm  (shared across all depths)\n', d1_fit);
    fprintf('RMSE = %.1f cm s^{-2}  |  R2 = %.4f\n', rmse, r2);

    for zi = 1:nZ
        fprintf('  z=%.2f cm:  F/m = %.1f ± %.1f cm s^{-2}\n', ...
            z_targets(zi), Fz_over_m(zi), Fz_over_m_se(zi));
    end

    fprintf('-------------------------------------------------------------\n\n');

    %% ── Bin display data (v >= v_min only) ───────────────────────────────
    nonempty = ~cellfun(@isempty, v_at_z);

    if ~any(nonempty)
        error('No nonempty depth groups after cutoff.');
    end

    v_global_max = max(cellfun(@max, v_at_z(nonempty)));
    v_edges      = linspace(v_min, v_global_max*1.02, nBins+1);
    v_centers    = (v_edges(1:end-1) + v_edges(2:end)) / 2;

    bin_mean = nan(nZ, nBins);
    bin_std  = nan(nZ, nBins);

    for zi = 1:nZ

        if isempty(v_at_z{zi})
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

    %% ── Figure ───────────────────────────────────────────────────────────
    fig = figure('Name','a+g vs v (Katsuragi)', ...
                 'ToolBar','none','MenuBar','none');

    fig.Position = [100 100 740 560];

    ax = axes(fig, 'Position', [0.11 0.12 0.82 0.83]);
    hold(ax,'on');

    % Fit lines start from 0 across full range
    v_line = linspace(0, v_global_max*1.04, 400);

    for zi = 1:nZ

        if isempty(v_at_z{zi})
            continue;
        end

        col = depth_cols(zi,:);

        % Faint raw scatter
        scatter(ax, v_at_z{zi}, ag_at_z{zi}, 5, col, 'filled', ...
            'MarkerFaceAlpha', 0.10, ...
            'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');

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

        % Shared-d1 fit line
        ag_fit_line = Fz_over_m(zi) + v_line.^2 ./ d1_fit;

        plot(ax, v_line, ag_fit_line, 'k--', ...
            'LineWidth', 1.8, ...
            'HandleVisibility', 'off');
    end

    tick_start = ceil(v_min / 50) * 50;

    set(ax, ...
        'FontSize',14, ...
        'Box','on', ...
        'LineWidth',1.2, ...
        'XColor',[0 0 0], ...
        'YColor',[0 0 0], ...
        'XMinorTick','on', ...
        'YMinorTick','on', ...
        'TickDir','in', ...
        'XLim',  [v_min * 0.9, v_global_max * 1.05], ...
        'XTick', tick_start : 50 : ceil(v_global_max/50)*50, ...
        'YLim',  [0, ag_max_plot]);

    grid(ax,'off');

    xlabel(ax,'$v$  (cm s$^{-1}$)', ...
        'FontSize',20, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);

    ylabel(ax,'$a+g$  (cm s$^{-2}$)', ...
        'FontSize',20, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);

    lgd = legend(ax,'show', ...
        'FontSize',11, ...
        'Box','on', ...
        'Interpreter','latex', ...
        'EdgeColor',[0.25 0.25 0.25], ...
        'Color',[1 1 1]);

    lgd.Location = 'southeast';

    text(ax, 0.97, 0.05, ...
        sprintf('$d_1 = %.2f$ cm\n$R^2 = %.3f$', d1_fit, r2), ...
        'FontSize',11, ...
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

    %% ── Inset: a+g vs v² — parallel lines prove shared d1 ────────────────
    ax_in = axes(fig, 'Position', [0.18 0.66 0.25 0.25]);
    hold(ax_in,'on');

    v2_max_all = max(cellfun(@(v) max(v.^2), v_at_z(nonempty)));
    v2_edges   = linspace(v_min^2, v2_max_all*1.02, nBins+1);
    v2_centers = (v2_edges(1:end-1) + v2_edges(2:end)) / 2;

    for zi = 1:nZ

        if isempty(v_at_z{zi})
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

        % Fit lines from 0
        v2_line     = linspace(0, v2_max_all*1.04, 300);
        ag_fit_line = Fz_over_m(zi) + v2_line ./ d1_fit;

        plot(ax_in, v2_line, ag_fit_line, 'k--', ...
            'LineWidth', 1.4, ...
            'HandleVisibility','off');
    end

    set(ax_in, ...
        'FontSize',10, ...
        'Box','on', ...
        'LineWidth',1.0, ...
        'XColor',[0 0 0], ...
        'YColor',[0 0 0], ...
        'TickDir','in', ...
        'XMinorTick','on', ...
        'YMinorTick','on', ...
        'XLim',[0, v2_max_all*1.04], ...
        'YLim',[0, max(bin_mean(:),[],'omitnan')*1.15]);

    grid(ax_in,'off');

    xlabel(ax_in,'$v^2$  (cm$^2$ s$^{-2}$)', ...
        'FontSize',20, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);

    ylabel(ax_in,'$a+g$  (cm s$^{-2}$)', ...
        'FontSize',20, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);
end