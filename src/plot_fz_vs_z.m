function [fig, k_over_m, fitStats] = plot_fz_vs_z(heights, d1, cmap, Fz_over_m, z_targets, d0)
% PLOT_FZ_VS_Z  Katsuragi & Durian (2007) Fig. 3b — exact approach.
%
%   Shows:
%       1. Colored curves: smoothed F(z)/m = (a+g) - v²/d1 per height group
%       2. Open symbols:   F(z_i)/m intercepts from plot_ag_vs_v2
%       3. Dotted line:    forced-origin fit  F/m = (k/m)·z
%                          (fitted to stable depths only, z >= z_min_fit)
%       4. Dashed line:    Tsimring-Volfson candidate (if d0 supplied)
%
%   Inputs:
%       heights    - struct array from group_trials_by_height
%       d1         - shared inertial length scale from plot_ag_vs_v2 [cm]
%       cmap       - [nH x 3] colormap matching heights order
%       Fz_over_m  - depth-specific intercepts from plot_ag_vs_v2 [cm s^-2]
%       z_targets  - depths corresponding to Fz_over_m [cm]
%       d0         - (optional) depth scale for Tsimring-Volfson curve [cm]
%
%   Outputs:
%       fig        - figure handle
%       k_over_m   - forced-origin linear slope [s^-2]
%       fitStats   - fit diagnostics struct

    if nargin < 6, d0 = []; end

    nH        = numel(heights);
    z_targets = z_targets(:);
    Fz_over_m = Fz_over_m(:);
    tol       = 0.05;    % depth window from plot_ag_vs_v2 [cm]
    z_min_fit = 1.2;     % exclude geometric transition zone from k/m fit
    sg_span   = 9;       % Savitzky-Golay window for F(z)/m curve smoothing
    clip_frac = 0.90;    % clip each curve at 90% of its max depth (avoid stop artifacts)

    % ── Valid intercept points ─────────────────────────────────────────────
    ok    = isfinite(z_targets) & isfinite(Fz_over_m) & z_targets > 0 & Fz_over_m > 0;
    z_pts = z_targets(ok);
    F_pts = Fz_over_m(ok);

    if numel(z_pts) < 2
        error('Need at least two valid F(z_i)/m points.');
    end

    % ── Forced-origin fit — stable depths only (z >= z_min_fit) ───────────
    fit_ok = z_pts >= z_min_fit;
    z_fit  = z_pts(fit_ok);
    F_fit  = F_pts(fit_ok);
    if numel(z_fit) < 2
        warning('Fewer than 2 stable-depth points; using all depths.');
        z_fit = z_pts; F_fit = F_pts;
    end

    k_over_m = z_fit \ F_fit;
    F_pred   = k_over_m .* z_fit;
    resid    = F_fit - F_pred;
    rmse     = sqrt(mean(resid.^2));
    ss_res   = sum(resid.^2);
    ss_tot   = sum((F_fit - mean(F_fit)).^2);
    r2       = 1 - ss_res / ss_tot;

    fitStats.k_over_m  = k_over_m;
    fitStats.rmse      = rmse;
    fitStats.r2        = r2;
    fitStats.z_pts     = z_pts;
    fitStats.F_pts     = F_pts;
    fitStats.z_fit     = z_fit;
    fitStats.F_fit     = F_fit;
    fitStats.z_min_fit = z_min_fit;
    fitStats.d1        = d1;
    fitStats.d0        = d0;

    fprintf('\n-- F(z)/m vs z (Fig. 3b) ------------------------------------\n');
    fprintf('d1        = %.3f cm\n', d1);
    fprintf('z_min_fit = %.2f cm  (geometric transition excluded)\n', z_min_fit);
    fprintf('k/m       = %.1f s^{-2}  (forced-origin, z >= %.2f cm)\n', k_over_m, z_min_fit);
    fprintf('RMSE      = %.1f cm s^{-2}\n', rmse);
    fprintf('R^2       = %.4f\n', r2);
    if ~isempty(d0)
        fprintf('d0        = %.3f cm  (Tsimring-Volfson)\n', d0);
    end
    fprintf('-------------------------------------------------------------\n\n');

    % ── Compute smoothed F(z)/m trajectory per height group ───────────────
    z_max_data = max(z_pts) * 1.15;
    z_ref      = linspace(0, z_max_data, 300);
    Fz_curves  = nan(nH, numel(z_ref));

    for j = 1:nH
        Fz_mat = [];
        for i = 1:heights(j).nTrials
            k   = heights(j).trials(i).kinematics;
            idx = k.impact_index:k.stopFrame;
            z_tr  = k.z_smooth(idx);
            ag_tr = k.a_plus_g(idx);
            v_tr  = k.v_smooth(idx);
            fz_tr = ag_tr - v_tr.^2 ./ d1;

            % Clip at clip_frac of max depth to avoid stop-frame artifacts
            z_stop = max(z_tr, [], 'omitnan');
            valid  = isfinite(z_tr) & isfinite(fz_tr) & ...
                     z_tr >= 0 & fz_tr > 0 & z_tr <= clip_frac * z_stop;
            if sum(valid) < sg_span + 2, continue; end

            % SG smooth to reduce noise before interpolation
            fz_sm = sgolayfilt(fz_tr(valid), 2, sg_span);

            [z_u, ia] = unique(z_tr(valid));
            fz_u = fz_sm(ia);
            Fz_mat = [Fz_mat; interp1(z_u, fz_u, z_ref, 'pchip', NaN)];
        end
        if ~isempty(Fz_mat)
            Fz_curves(j,:) = mean(Fz_mat, 1, 'omitnan');
        end
    end

    % ── Figure ────────────────────────────────────────────────────────────
    fig = figure('Name','F(z)/m vs z  (Fig. 3b)', ...
                 'ToolBar','none','MenuBar','none');
    fig.Position = [100 100 660 520];
    ax = axes(fig, 'Position', [0.13 0.13 0.78 0.82]);
    hold(ax,'on');

    % 1. Colored curves — one per height group
    for j = 1:nH
        good = isfinite(Fz_curves(j,:));
        if ~any(good), continue; end
        plot(ax, z_ref(good), Fz_curves(j,good), '-', ...
            'Color', cmap(j,:), 'LineWidth', 1.8, ...
            'HandleVisibility', 'off');
    end

    % 2. Dotted line: k/m fitted to stable depths
    plot(ax, z_ref, k_over_m .* z_ref, 'k:', 'LineWidth', 2.4, ...
        'DisplayName', sprintf('$(k/m)z,\\;k/m=%.0f$ s$^{-2}$', k_over_m));

    % 3. Tsimring-Volfson dashed curve (optional)
    if ~isempty(d0)
        g_cm = 980;
        F_tv = g_cm .* (1 + (3*(z_ref./d0).^2 - 1) .* exp(-2.*z_ref./d1));
        plot(ax, z_ref, F_tv, '--', 'Color', [0.50 0.50 0.50], ...
            'LineWidth', 1.8, ...
            'DisplayName', '$g\{1+[3(z/d_0)^2-1]e^{-2z/d_1}\}$');
    end

    % 4. Open symbols: all intercepts (stable = filled open circle, transition = grey)
    for ii = 1:numel(z_pts)
        if z_pts(ii) >= z_min_fit
            mfc = 'w';   col_pt = 'k';    % stable — open black circle
        else
            mfc = 'w';   col_pt = [0.6 0.6 0.6];  % transition — open grey circle
        end
        errorbar(ax, z_pts(ii), F_pts(ii), 0, 0, tol, tol, ...
            'o', 'Color', col_pt, ...
            'MarkerFaceColor', mfc, 'MarkerEdgeColor', col_pt, ...
            'MarkerSize', 9, 'LineWidth', 1.6, 'CapSize', 4, ...
            'HandleVisibility', 'off');
    end
    % Dummy entry for legend
    plot(ax, nan, nan, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineWidth', 1.6, ...
        'DisplayName', '$F(z_i)/m$ intercepts');

    % ── Axes ──────────────────────────────────────────────────────────────
    y_max = max([F_pts; Fz_curves(:)], [], 'omitnan') * 1.15;
    if ~isfinite(y_max) || y_max <= 0, y_max = max(F_pts) * 1.25; end

    set(ax, 'FontSize',13, 'Box','on', 'LineWidth',1.2, ...
        'XColor',[0 0 0], 'YColor',[0 0 0], ...
        'XMinorTick','on', 'YMinorTick','on', 'TickDir','in', ...
        'XLim',[0, z_max_data], 'YLim',[0, y_max]);
    grid(ax,'off');

    xlabel(ax, '$z$  (cm)', 'FontSize',16, 'Interpreter','latex', 'Color',[0 0 0]);
    ylabel(ax, '$F(z)/m$  (cm s$^{-2}$)', 'FontSize',16, 'Interpreter','latex', 'Color',[0 0 0]);

    legend(ax, 'show', 'FontSize',10, 'Interpreter','latex', ...
        'Box','on', 'EdgeColor',[0.25 0.25 0.25], 'Location','northwest');

    text(ax, 0.97, 0.05, ...
        sprintf('$d_1 = %.2f$ cm\n$k/m = %.0f$ s$^{-2}$\n$R^2 = %.3f$\n$(z \\geq %.1f$ cm$)$', ...
            d1, k_over_m, r2, z_min_fit), ...
        'Units','normalized', 'Interpreter','latex', 'FontSize',10.5, ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
        'BackgroundColor',[1 1 1], 'EdgeColor',[0.25 0.25 0.25], 'Margin',4);
end