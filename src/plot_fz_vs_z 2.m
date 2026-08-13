function [fig, k_over_m, fitStats] = plot_fz_vs_z(heights, d1, ~, ...
        z_targets, Fz_over_m, Fz_per_height, Fz_per_height_se, ...
        v0_per_height, d0)
% PLOT_FZ_VS_Z  Step 1 — F(z)/m vs z  (Katsuragi Fig. 3b comparison).
%
%   Symbol style:
%     Black edge, hollow colored fill by regime:
%       near-surface (z < 0.40 cm) : light grey fill
%       transition   (0.40-1.30 cm): muted mauve fill
%       stable shaft (z >= 1.30 cm): dark crimson-pink fill
%     Vertical bars: OLS-based SE at all depths
%     Horizontal bars: +/-0.05 cm (OLS window half-width)
%
%   Reference line colors (consistent across Step figures):
%     k|z|        : black dotted
%     F~V         : dark blue dashed  (plot_fz_vs_V)
%     F~A(z)z     : dark green solid  (plot_fz_vs_Az)
%     F~[A(z)z]^p : dark red dashed   (plot_fz_vs_Az)

    if nargin < 9, d0 = []; end

    z_targets    = z_targets(:);
    Fz_over_m    = Fz_over_m(:);
    nZ           = numel(z_targets);
    z_min_fit    = 1.30;
    z_detect     = 0.40;
    tol_disp     = 0.05;   % horizontal error bar half-width [cm]

    % ── Regime fill colors (black edge on all) ────────────────────────────
    col_near   = [0.82 0.82 0.82];   % light grey
    col_trans  = [0.55 0.35 0.45];   % muted mauve
    col_stable = [0.55 0.05 0.25];   % dark crimson-pink
    ms = 9; lw_sym = 1.6;

    %% ── Compute OLS-based SE at all depths ───────────────────────────────
    % SE of OLS intercept b0:
    %   SE_b0 = sqrt(MSE * (1/n + v2_mean^2 / S_v2v2))
    % where S_v2v2 = sum((v2 - v2_mean)^2), MSE = RSS/(n-2)
    % Requires re-collecting v data per depth.

    Fz_ols_se = nan(nZ, 1);
    tol_collect = 0.05;
    ag_cap      = 60000;
    v_min       = 40;

    for zi = 1:nZ
        v_zi = []; ag_zi = [];
        for j = 1:numel(heights)
            for i = 1:heights(j).nTrials
                k_t  = heights(j).trials(i).kinematics;
                idx  = k_t.impact_index:k_t.stopFrame;
                z_t  = k_t.z_smooth(idx);
                v_t  = k_t.v_smooth(idx);
                ag_t = k_t.a_plus_g(idx);
                in   = abs(z_t - z_targets(zi)) < tol_collect & ...
                       v_t >= v_min & isfinite(v_t) & isfinite(ag_t) & ...
                       v_t > 0 & ag_t > 0 & ag_t < ag_cap;
                v_zi  = [v_zi;  v_t(in)];
                ag_zi = [ag_zi; ag_t(in)];
            end
        end
        n_zi = numel(v_zi);
        if n_zi < 3, continue; end
        v2_zi    = v_zi.^2;
        A_zi     = [ones(n_zi,1), v2_zi];
        coeffs   = A_zi \ ag_zi;
        ag_pred  = A_zi * coeffs;
        RSS      = sum((ag_zi - ag_pred).^2);
        MSE      = RSS / (n_zi - 2);
        v2_mean  = mean(v2_zi);
        Svv      = sum((v2_zi - v2_mean).^2);
        if Svv > 0
            Fz_ols_se(zi) = sqrt(MSE * (1/n_zi + v2_mean^2/Svv));
        end
    end

    %% ── SE for stable depths from per-height spread ──────────────────────
    % Use per-height spread where available (more honest unit)
    Fz_pooled_se = nan(nZ, 1);
    for zi = 1:nZ
        if isempty(Fz_per_height), break; end
        vals = Fz_per_height(:, zi);
        ok   = isfinite(vals);
        if nnz(ok) >= 2
            Fz_pooled_se(zi) = std(vals(ok)) / sqrt(nnz(ok));
        end
    end

    % Use per-height SE for stable depths, OLS SE elsewhere
    SE_plot = nan(nZ, 1);
    for zi = 1:nZ
        if z_targets(zi) >= z_min_fit && isfinite(Fz_pooled_se(zi))
            SE_plot(zi) = Fz_pooled_se(zi);
        elseif isfinite(Fz_ols_se(zi))
            SE_plot(zi) = Fz_ols_se(zi);
        else
            SE_plot(zi) = 0;
        end
    end

    %% ── k/m fit from stable depths ───────────────────────────────────────
    is_stable = z_targets >= z_min_fit & isfinite(Fz_over_m) & Fz_over_m > 0;
    z_stable  = z_targets(is_stable);
    F_stable  = Fz_over_m(is_stable);
    SE_stable = SE_plot(is_stable);

    z_fit = []; F_fit = [];
    if ~isempty(Fz_per_height)
        for zi = find(is_stable)'
            vals = Fz_per_height(:, zi);
            ok_j = isfinite(vals) & vals > 0;
            z_fit = [z_fit; repmat(z_targets(zi), nnz(ok_j), 1)];
            F_fit = [F_fit; vals(ok_j)];
        end
    end
    if numel(z_fit) < 2, z_fit = z_stable; F_fit = F_stable; end

    k_over_m = z_fit \ F_fit;
    resid    = F_fit - k_over_m .* z_fit;
    r2_lin   = 1 - sum(resid.^2) / sum((F_fit - mean(F_fit)).^2);
    rmse_lin = sqrt(mean(resid.^2));

    fprintf('\n-- plot_fz_vs_z: k|z| fit (z >= %.2f cm) -------------------\n', z_min_fit);
    fprintf('k/m  = %.1f s^-2  |  R^2 = %.4f\n', k_over_m, r2_lin);
    fprintf('-------------------------------------------------------------\n\n');

    %% ── Axis limits ──────────────────────────────────────────────────────
    all_F_valid = Fz_over_m(isfinite(Fz_over_m));
    y_lo = min(0, min(all_F_valid) * 1.15);
    y_hi = max(F_stable + SE_stable, [], 'omitnan') * 1.25;
    if ~isfinite(y_hi) || y_hi <= 0, y_hi = max(F_stable)*1.25; end
    z_max_data = max(z_targets) * 1.12;

    %% ── Figure ───────────────────────────────────────────────────────────
    fig = figure('Name','Step 1: F(z)/m vs z', ...
                 'ToolBar','none','MenuBar','none');
    fig.Position = [100 100 760 600];
    ax = axes(fig, 'Position', [0.12 0.12 0.85 0.83]);
    hold(ax, 'on');

    % y=0 reference — black dashed
    yline(ax, 0, '--', 'Color',[0.10 0.10 0.10], 'LineWidth',1.8, ...
        'HandleVisibility','off');

    % k|z| dotted — black (Step 1 reference color)
    z_line = linspace(0, z_max_data, 300);
    plot(ax, z_line, k_over_m .* z_line, 'k:', 'LineWidth', 3.0, ...
        'DisplayName', sprintf('(k/m)z,  k/m = %.0f s^{-2}', k_over_m));

    % Transition boundary — dark red dashed vertical
    xline(ax, z_min_fit, '--', 'Color',[0.65 0.05 0.05], 'LineWidth', 2.5, ...
        'HandleVisibility','off');
    text(ax, z_min_fit - 0.04, y_hi * 0.88, ...
        sprintf('geometric\ntransition'), ...
        'FontSize', 18, 'Interpreter', 'none', ...
        'Color', [0.65 0.05 0.05], ...
        'HorizontalAlignment','right', 'VerticalAlignment','top');
    text(ax, z_min_fit + 0.04, y_hi * 0.88, ...
        sprintf('stable shaft\nz >= %.1f cm', z_min_fit), ...
        'FontSize', 18, 'Interpreter', 'none', ...
        'Color', [0.65 0.05 0.05], ...
        'HorizontalAlignment','left', 'VerticalAlignment','top');

    % Plot all symbols
    for zi = 1:nZ
        if ~isfinite(Fz_over_m(zi)), continue; end
        z_q  = z_targets(zi);
        F_q  = Fz_over_m(zi);
        SE_q = SE_plot(zi);

        if z_q < z_detect
            fc = col_near;
        elseif z_q < z_min_fit
            fc = col_trans;
        else
            fc = col_stable;
        end

        errorbar(ax, z_q, F_q, SE_q, SE_q, tol_disp, tol_disp, ...
            'o', ...
            'Color',           'k', ...
            'MarkerFaceColor', fc, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize',      ms, ...
            'LineWidth',       lw_sym, ...
            'CapSize',         3, ...
            'HandleVisibility','off');
    end


    set(ax, 'FontSize',13, 'Box','on', 'LineWidth',1.2, ...
        'XColor',[0 0 0], 'YColor',[0 0 0], ...
        'XMinorTick','on', 'YMinorTick','on', 'TickDir','in', ...
        'XLim',[0, z_max_data], 'YLim',[y_lo, y_hi]);
    grid(ax,'off');
    xlabel(ax, '$z$  (cm)', 'FontSize',16, 'Interpreter','latex', 'Color',[0 0 0]);
    ylabel(ax, '$F(z)/m$  (cm s$^{-2}$)', 'FontSize',16, 'Interpreter','latex', 'Color',[0 0 0]);


    %% ── Pack fitStats ────────────────────────────────────────────────────
    fitStats.k_over_m         = k_over_m;
    fitStats.r2_lin           = r2_lin;
    fitStats.rmse_lin         = rmse_lin;
    fitStats.z_min_fit        = z_min_fit;
    fitStats.z_detect         = z_detect;
    fitStats.d1               = d1;
    fitStats.z_stable         = z_stable;
    fitStats.F_stable         = F_stable;
    fitStats.SE_stable        = SE_stable;
    fitStats.Fz_over_m        = Fz_over_m;
    fitStats.Fz_ols_se        = Fz_ols_se;
    fitStats.SE_plot          = SE_plot;
    fitStats.Fz_pooled_se     = Fz_pooled_se;
    fitStats.Fz_per_height    = Fz_per_height;
    fitStats.Fz_per_height_se = Fz_per_height_se;
    fitStats.v0_per_height    = v0_per_height;
end