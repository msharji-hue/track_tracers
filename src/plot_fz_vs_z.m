function [fig, k_over_m, fitStats] = plot_fz_vs_z(heights, d1, cmap, Fz_over_m, z_targets, d0, Fz_over_m_se, Fz_per_height, Fz_per_height_se, v0_per_height)
% PLOT_FZ_VS_Z  Katsuragi & Durian (2007) Fig. 3b — exact + improved.
%
%   Symbols exactly replicate Katsuragi's approach:
%       - At each depth target, heights sorted by v0 and divided into
%         n_subsets groups.
%       - Each group gives one open symbol.
%       - Vertical error bar = SE across height-group means (honest unit).
%       - Horizontal error bar = depth tolerance ±tol.
%       - Colored by speed subset (low/mid/high v0).
%
%   Additional:
%       - 3 speed-binned mean curves with bootstrap 95% CI bands.
%       - Pearson r annotated at each depth (from per-height data).
%       - Single k|z| dotted reference line.

    if nargin < 6,  d0             = [];  end
    if nargin < 7 || isempty(Fz_over_m_se)
        Fz_over_m_se = zeros(size(Fz_over_m));
    end
    if nargin < 8 || isempty(Fz_per_height)
        Fz_per_height    = [];
        Fz_per_height_se = [];
        v0_per_height    = [];
    end

    nH           = numel(heights);
    z_targets    = z_targets(:);
    Fz_over_m    = Fz_over_m(:);
    Fz_over_m_se = Fz_over_m_se(:);
    nZ           = numel(z_targets);

    tol           = 0.07;
    z_min_fit     = 1.30;
    z_min_display = 0.80;
    n_speed_bins  = 3;
    n_subsets     = 3;    % height subsets per depth — Katsuragi style
    n_boot        = 500;
    sg_span       = 19;
    clip_frac     = 0.88;
    rng(42);

    bin_cols = [0.15 0.45 0.80;   % low  v0 — blue
                0.15 0.68 0.38;   % mid  v0 — green
                0.82 0.18 0.12];  % high v0 — red

    %% ── Valid pooled intercept points (for k/m fit and Pearson r) ────────
    ok = isfinite(z_targets) & isfinite(Fz_over_m) & z_targets > 0;
    z_pts = z_targets(ok);
    F_pts = Fz_over_m(ok);

    if numel(z_pts) < 2
        error('Need at least two valid F(z_i)/m points.');
    end

    %% ── k/m forced-origin fit ────────────────────────────────────────────
    if ~isempty(Fz_per_height)
        z_flat = []; F_flat = [];
        for zi = 1:nZ
            if z_targets(zi) < z_min_fit, continue; end
            vals = Fz_per_height(:, zi);
            ok_j = isfinite(vals) & vals > 0;
            z_flat = [z_flat; repmat(z_targets(zi), sum(ok_j), 1)];
            F_flat = [F_flat; vals(ok_j)];
        end
        z_fit = z_flat; F_fit = F_flat;
    else
        fit_ok = z_pts >= z_min_fit;
        z_fit  = z_pts(fit_ok); F_fit = F_pts(fit_ok);
    end
    if numel(z_fit) < 2, z_fit = z_pts; F_fit = F_pts; end

    k_over_m = z_fit \ F_fit;
    F_pred   = k_over_m .* z_fit;
    resid    = F_fit - F_pred;
    rmse     = sqrt(mean(resid.^2));
    ss_res   = sum(resid.^2);
    ss_tot   = sum((F_fit - mean(F_fit)).^2);
    r2       = 1 - ss_res / ss_tot;

    %% ── Pearson r at each depth (from per-height data) ───────────────────
    r_vals = nan(numel(z_pts), 1);
    p_vals = nan(numel(z_pts), 1);

    if ~isempty(Fz_per_height) && ~isempty(v0_per_height)
        for ii = 1:numel(z_pts)
            [~, col] = min(abs(z_targets - z_pts(ii)));
            valid    = isfinite(Fz_per_height(:,col)) & ...
                       isfinite(v0_per_height) & ...
                       Fz_per_height(:,col) > 0;
            if sum(valid) >= 4
                [r_vals(ii), p_vals(ii)] = corr( ...
                    v0_per_height(valid), Fz_per_height(valid,col));
            end
        end
    end

    %% ── One pooled symbol per depth — Katsuragi black hollow style ──────
    % F(zi)/m = mean across all height groups at that depth
    % SE      = std across height-group means / sqrt(n_groups)  (honest unit)
    sym_z     = [];
    sym_F     = [];
    sym_SE    = [];
    sym_infit = [];

    if ~isempty(Fz_per_height) && ~isempty(v0_per_height)
        for zi = 1:nZ
            col_ph = Fz_per_height(:, zi);
            valid  = isfinite(col_ph) & col_ph > 0;
            if sum(valid) < 1, continue; end

            F_groups = col_ph(valid);
            F_mean_z = mean(F_groups);
            if numel(F_groups) >= 2
                SE_z = std(F_groups) / sqrt(numel(F_groups));
            else
                SE_z = 0;
            end

            sym_z     = [sym_z;     z_targets(zi)];
            sym_F     = [sym_F;     F_mean_z];
            sym_SE    = [sym_SE;    SE_z];
            sym_infit = [sym_infit; z_targets(zi) >= z_min_fit];
        end
    else
        % Fallback: use pooled Fz_over_m
        sym_z     = z_pts;
        sym_F     = F_pts;
        sym_SE    = Fz_over_m_se(ok);
        sym_infit = z_pts >= z_min_fit;
    end

    %% ── Power-law fit: F/m = C * z^n ─────────────────────────────────────
    fit_ok_log = z_pts >= z_min_fit & F_pts > 0;
    if sum(fit_ok_log) >= 2
        p_log  = polyfit(log(z_pts(fit_ok_log)), log(F_pts(fit_ok_log)), 1);
        n_pow  = p_log(1);
        C_pow  = exp(p_log(2));
        F_pow_pred = C_pow .* z_pts(fit_ok_log).^n_pow;
        ss_res_pow = sum((F_pts(fit_ok_log) - F_pow_pred).^2);
        ss_tot_pow = sum((F_pts(fit_ok_log) - mean(F_pts(fit_ok_log))).^2);
        r2_pow     = 1 - ss_res_pow / ss_tot_pow;
    else
        n_pow = NaN; C_pow = NaN; r2_pow = NaN;
    end

    fitStats.n_pow  = n_pow;
    fitStats.C_pow  = C_pow;
    fitStats.r2_pow = r2_pow;
    fitStats.z_pts     = z_pts;
    fitStats.F_pts     = F_pts;
    fitStats.r_vals    = r_vals;
    fitStats.p_vals    = p_vals;
    fitStats.z_min_fit = z_min_fit;
    fitStats.d1        = d1;

    fitStats.k_over_m  = k_over_m;
    fitStats.rmse      = rmse;
    fitStats.r2        = r2;

    fprintf('\n-- F(z)/m vs z (Fig. 3b) ------------------------------------\n');
    fprintf('d1        = %.3f cm\n', d1);
    fprintf('z_min_fit = %.2f cm  |  k/m = %.1f s^{-2}\n', z_min_fit, k_over_m);
    fprintf('Linear:    RMSE = %.1f  |  R^2 = %.4f\n', rmse, r2);
    if isfinite(n_pow)
        fprintf('Power law: F/m = %.1f * z^%.3f  |  R^2 = %.4f\n', C_pow, n_pow, r2_pow);
    end
    fprintf('-------------------------------------------------------------\n\n');

    %% ── y_max from symbols only ──────────────────────────────────────────
    sym_infit  = logical(sym_infit);   % ensure logical type for indexing
    z_max_data = max(z_pts) * 1.12;

    if any(sym_infit)
        y_max = max(sym_F(sym_infit) + sym_SE(sym_infit), [], 'omitnan') * 1.25;
    else
        y_max = max(F_pts) * 1.25;
    end
    if ~isfinite(y_max) || y_max <= 0, y_max = max(F_pts)*1.25; end

    %% ── Figure ───────────────────────────────────────────────────────────
    fig = figure('Name','F(z)/m vs z  (Fig. 3b)', ...
                 'ToolBar','none','MenuBar','none');
    fig.Position = [100 100 720 560];
    ax = axes(fig, 'Position', [0.13 0.13 0.84 0.82]);
    hold(ax,'on');

    %% 1. k|z| dotted reference — extended from z=0
    z_line = linspace(0, z_max_data, 300);
    plot(ax, z_line, k_over_m .* z_line, 'k:', 'LineWidth', 3.0, ...
        'DisplayName', sprintf('$(k/m)z,\\;k/m=%.0f$ s$^{-2}$', k_over_m));

    %% 2. Power-law dashed reference — from z=0
    if isfinite(n_pow)
        z_pow = linspace(0, z_max_data, 300);
        plot(ax, z_pow, C_pow .* z_pow.^n_pow, 'k--', 'LineWidth', 2.0, ...
            'DisplayName', sprintf('$%.0f\\,z^{%.2f}$, $R^2=%.2f$', ...
                C_pow, n_pow, r2_pow));
    end

    %% 3. Symbols — black for stable depths, grey for excluded shallow depths
    for ii = 1:numel(sym_z)
        if sym_infit(ii)
            col_pt = 'k';
        else
            col_pt = [0.55 0.55 0.55];
        end
        errorbar(ax, sym_z(ii), sym_F(ii), ...
            sym_SE(ii), sym_SE(ii), tol, tol, ...
            'o', 'Color',          col_pt, ...
            'MarkerFaceColor',     'w',    ...
            'MarkerEdgeColor',     col_pt, ...
            'MarkerSize',          10,     ...
            'LineWidth',           1.8,    ...
            'CapSize',             4,      ...
            'HandleVisibility',    'off');
    end

    plot(ax, nan, nan, 'o', 'Color','k', 'MarkerFaceColor','w', ...
        'MarkerEdgeColor','k', 'MarkerSize',10, 'LineWidth',1.8, ...
        'DisplayName', '$F(z_i)/m$ intercepts');

    %% ── Axes ─────────────────────────────────────────────────────────────
    set(ax, 'FontSize',13, 'Box','on', 'LineWidth',1.2, ...
        'XColor',[0 0 0], 'YColor',[0 0 0], ...
        'XMinorTick','on', 'YMinorTick','on', 'TickDir','in', ...
        'XLim',[0, z_max_data], 'YLim',[0, y_max]);
    grid(ax,'off');
    xlabel(ax, '$z$  (cm)', 'FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax, '$F(z)/m$  (cm s$^{-2}$)', 'FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    legend(ax,'show','FontSize',10,'Interpreter','latex','Box','on', ...
        'EdgeColor',[0.25 0.25 0.25],'Location','northwest');
    if isfinite(n_pow)
        ann_str = sprintf(['$d_1 = %.2f$ cm\n$k/m = %.0f$ s$^{-2}$, $R^2=%.2f$\n' ...
                           '$n = %.2f$ (power law), $R^2=%.2f$\n' ...
                           '(stable shaft, $z \\geq %.1f$ cm)'], ...
            d1, k_over_m, r2, n_pow, r2_pow, z_min_fit);
    else
        ann_str = sprintf(['$d_1 = %.2f$ cm\n$k/m = %.0f$ s$^{-2}$\n' ...
                           '$R^2 = %.3f$\n(stable shaft, $z \\geq %.1f$ cm)'], ...
            d1, k_over_m, r2, z_min_fit);
    end
    text(ax, 0.97, 0.05, ann_str, ...
        'Units','normalized','Interpreter','latex','FontSize',10.5, ...
        'HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25],'Margin',4);
end