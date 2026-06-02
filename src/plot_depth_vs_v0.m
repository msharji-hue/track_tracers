function [fig, d0] = plot_depth_vs_v0(heights, cmap, ~, d1, k_over_m)
% PLOT_DEPTH_VS_V0  Final penetration depth vs impact speed.
%   Replicates Katsuragi & Durian 2007 Fig. 2b.
%
%   Always fits and returns d0 from the (d0²H)^(1/3) depth-scaling law.
%
%   Forward model overlay (Katsuragi dotted curve) shown only when
%   both d1 and k_over_m are supplied. Uses run_forward_model_simple —
%   no area scaling, direct analogue of their approach.
%
%   Inputs:
%       heights   - struct array from group_trials_by_height
%       cmap      - [nH x 3] colormap
%       ~         - ignored (footArea no longer used)
%       d1        - kinematic inertial length scale [cm] (from Fig 3a)
%       k_over_m  - friction coefficient [s^-2] (from Fig 3b)
%
%   Outputs:
%       fig  - figure handle
%       d0   - depth scale from (d0²H)^(1/3) fit [cm] — feeds plot_fz_vs_z

    doForward = nargin >= 4 && ~isempty(d1) && ...
                nargin >= 5 && ~isempty(k_over_m);

    g_cm = 980;

    % ── Collect per-trial scalars ─────────────────────────────────────────
    v0_all = [];
    d_all  = [];
    H_all  = [];

    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials
            s  = heights(j).trials(i).scalars;
            h  = heights(j).h_cm;
            d  = s.d_final_cm;
            v0_all(end+1) = s.v0_cm_s;
            d_all(end+1)  = d;
            H_all(end+1)  = h + d;
        end
    end

    ss_tot = sum((d_all - mean(d_all)).^2);

    % ── Linear fit: d = d0_lin + alpha*|v0| ──────────────────────────────
    X         = [ones(numel(v0_all),1), v0_all(:)];
    coeffs    = X \ d_all(:);
    d0_lin    = coeffs(1);
    alpha_fit = coeffs(2);
    r2_lin    = 1 - sum((d_all - (d0_lin + alpha_fit.*v0_all)).^2) / ss_tot;

    % ── Depth-scaling fit: d = (d0² * H)^(1/3) ───────────────────────────
    d0_sq = mean(d_all.^3 ./ H_all);
    d0    = sqrt(d0_sq);

    d_scaling_pred = (d0^2 .* H_all).^(1/3);
    r2_scaling     = 1 - sum((d_all - d_scaling_pred).^2) / ss_tot;

    fprintf('\n-- d vs v0 --------------------------------------------------\n');
    fprintf('Linear:         d0=%.3f cm  alpha=%.5f s/cm  R2=%.4f\n', ...
        d0_lin, alpha_fit, r2_lin);
    fprintf('Depth-scaling:  d0=%.3f cm  (d0^2*H)^1/3     R2=%.4f\n', ...
        d0, r2_scaling);

    % ── Forward model (optional) ──────────────────────────────────────────
    v0_line = linspace(min(v0_all)*0.90, max(v0_all)*1.05, 80);
    d_fwd   = nan(size(v0_line));
    r2_fwd  = NaN;

    if doForward
        % Vectorized integration — all v0 values stepped simultaneously.
        % ~100x faster than calling run_forward_model_simple in a loop.
        dz        = 0.005;                    % step size [cm]
        max_depth = max(d_all) * 2.0;         % safety ceiling
        g         = g_cm;

        % Combine curve + trial v0 into one batch
        v0_batch  = [v0_line(:); v0_all(:)];
        nB        = numel(v0_batch);
        v2        = v0_batch.^2;
        d_batch   = zeros(nB, 1);
        stopped   = false(nB, 1);
        z_now     = 0;

        while ~all(stopped) && z_now < max_depth
            F        = g - k_over_m .* z_now - v2 ./ d1;
            v2_new   = v2 + 2 .* F .* dz;
            hit      = ~stopped & v2_new <= 0;
            d_batch(hit)  = z_now;
            stopped(hit)  = true;
            v2        = max(v2_new, 0);
            v2(stopped) = 0;
            z_now     = z_now + dz;
        end

        % Remaining (never stopped = hit max_depth ceiling)
        d_batch(~stopped) = max_depth;

        nL            = numel(v0_line);
        d_fwd         = d_batch(1:nL)';
        d_pred_trials = d_batch(nL+1:end)';

        r2_fwd = 1 - sum((d_all - d_pred_trials).^2) / ss_tot;
        fprintf('Forward model:  d1=%.3f cm  k/m=%.1f s^-2   R2=%.4f\n', ...
            d1, k_over_m, r2_fwd);
    end

    fprintf('-------------------------------------------------------------\n\n');

    % ── Reference curves ──────────────────────────────────────────────────
    d_lin_ref = d0_lin + alpha_fit .* v0_line;

    % ── Figure ────────────────────────────────────────────────────────────
    fig = figure('Name','d vs v0','ToolBar','none','MenuBar','none');
    fig.Position = [100 100 660 520];
    ax = axes(fig, 'Position', [0.12 0.13 0.60 0.82]);
    hold(ax,'on');

    % Linear reference — dashed grey (Katsuragi Fig. 2b grey dashed)
    plot(ax, v0_line, d_lin_ref, '--', 'Color', [0.25 0.25 0.25], ...
        'LineWidth', 2.0, 'HandleVisibility', 'off');

    % Forward model — dotted black (Katsuragi Fig. 2b dotted black)
    if doForward
        plot(ax, v0_line, d_fwd, 'k:', 'LineWidth', 2.0, ...
            'HandleVisibility', 'off');
    end

    % Data points per height group
    for j = 1:numel(heights)
        hg          = heights(j);
        [~, marker] = get_height_style(hg.h_cm);
        errorbar(ax, hg.v0_mean, hg.d_mean, hg.d_std, hg.d_std, ...
            hg.v0_std, hg.v0_std, marker, ...
            'Color', cmap(j,:), 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', cmap(j,:), 'MarkerSize', 9, 'LineWidth', 1.8, ...
            'HandleVisibility', 'off');
    end

    set(ax,'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in', ...
        'XLim',[min(v0_all)*0.90, max(v0_all)*1.05], ...
        'YLim',[0, max(d_all)*1.15]);
    grid(ax,'off');
    xlabel(ax,'$v_0$  (cm s$^{-1}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax,'$d$  (cm)',            'FontSize',16,'Interpreter','latex','Color',[0 0 0]);

    % ── Legend (dummy axes) ───────────────────────────────────────────────
    ax2 = axes(fig, 'Position', ax.Position, 'Visible','off');
    hold(ax2,'on');

    plot(ax2, nan, nan, '--', 'Color',[0.25 0.25 0.25], 'LineWidth',2.0, ...
        'DisplayName', sprintf('$d = %.2f + %.4f\\,v_0$', d0_lin, alpha_fit));

    if doForward
        plot(ax2, nan, nan, 'k:', 'LineWidth',2.0, ...
            'DisplayName', sprintf('Forward model ($d_1=%.2f$ cm, $k/m=%.0f$ s$^{-2}$), $R^2=%.3f$', ...
                d1, k_over_m, r2_fwd));
    end

    lgd2 = legend(ax2,'show', ...
        'FontSize',9, 'Box','on', 'Interpreter','latex', ...
        'EdgeColor',[0.15 0.15 0.15], 'Color',[1 1 1]);
    lgd2.Location = 'none';

    if doForward
        lgd2.Position = [0.735 0.60 0.245 0.14];
    else
        lgd2.Position = [0.57 0.16 0.34 0.08];
    end

    % ── Annotation ────────────────────────────────────────────────────────
    if doForward
        ann_str = sprintf(['$R^2_{\\rm lin} = %.3f$\n' ...
                           '$R^2_{\\rm fwd} = %.3f$\n' ...
                           '$d_0 = %.2f$ cm'], ...
            r2_lin, r2_fwd, d0);
    else
        ann_str = sprintf(['$R^2_{\\rm lin} = %.3f$\n' ...
                           '$d_0 = %.2f$ cm'], ...
            r2_lin, d0);
    end

    text(ax, 0.03, 0.97, ann_str, ...
        'FontSize',9,'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'Interpreter','latex','FontAngle','italic', ...
        'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25], ...
        'LineWidth',1.0,'Margin',4);
end
