function [fig, d0] = plot_depth_vs_v0(heights, cmap, ~, d1, k_over_m)
% PLOT_DEPTH_VS_V0  Final penetration depth vs impact speed.
%   Replicates Katsuragi & Durian 2007 Fig. 2b.
%
%   Three reference curves:
%     1. d = d0_lin + alpha*v0          (empirical linear, grey dashed)
%     2. d = (d0^2 * H(v0))^(1/3)      (depth-scaling, grey dash-dot)
%        H(v0) computed correctly as h(v0) + d_fwd(v0) — not mean H
%     3. Forward model ODE              (black dotted, requires d1+k_over_m)
%
%   Inputs:
%       heights   - struct array from group_trials_by_height
%       cmap      - [nH x 3] colormap
%       ~         - ignored
%       d1        - inertial length scale [cm] (from plot_ag_vs_v2)
%       k_over_m  - friction coefficient [s^-2] (from plot_fz_vs_z)
%
%   Outputs:
%       fig - figure handle
%       d0  - depth scale from (d0^2*H)^(1/3) fit [cm]

    doForward = nargin >= 4 && ~isempty(d1) && ...
                nargin >= 5 && ~isempty(k_over_m);

    g_cm = 980;

    %% ── Collect per-trial scalars ────────────────────────────────────────
    v0_all = [];
    d_all  = [];
    H_all  = [];

    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials
            s  = heights(j).trials(i).scalars;
            v0_all(end+1) = s.v0_cm_s;
            d_all(end+1)  = s.d_final_cm;
            H_all(end+1)  = heights(j).h_cm + s.d_final_cm;
        end
    end

    ss_tot = sum((d_all - mean(d_all)).^2);

    %% ── Depth-scaling fit: d = (d0^2*H)^(1/3) — fitted first ───────────
    % Katsuragi fits d0 from the depth-scaling law first, then uses the
    % SAME d0 as the intercept of the linear fit d = d0 + alpha*|v0|.
    % This ensures both curves share the same y-intercept at v0=0.
    d0_sq          = mean(d_all.^3 ./ H_all);
    d0             = sqrt(d0_sq);
    d_scaling_pred = (d0^2 .* H_all).^(1/3);
    r2_scaling     = 1 - sum((d_all-d_scaling_pred).^2)/ss_tot;

    %% ── Linear fit: d = d0 + alpha*v0 — d0 shared with depth-scaling ────
    % Force intercept = d0 from depth-scaling (Katsuragi Fig 2b method).
    % Only alpha is free — fit by minimizing sum((d - d0 - alpha*v0)^2).
    alpha_fit = v0_all(:) \ (d_all(:) - d0);
    r2_lin    = 1 - sum((d_all-(d0+alpha_fit.*v0_all)).^2)/ss_tot;
    d0_lin    = d0;   % shared intercept

    fprintf('\n-- d vs v0 --------------------------------------------------\n');
    fprintf('Depth-scaling:  d0=%.3f cm  (d0^2*H)^1/3     R2=%.4f\n', ...
        d0, r2_scaling);
    fprintf('Linear:         d0=%.3f cm  alpha=%.5f s/cm  R2=%.4f\n', ...
        d0_lin, alpha_fit, r2_lin);
    fprintf('  (d0 shared between both fits — matches Katsuragi Fig 2b)\n');

    %% ── Reference curves ─────────────────────────────────────────────────
    v0_max  = max(v0_all) * 1.08;
    v0_line = linspace(0, v0_max, 200);

    % 1. Linear reference
    d_lin_ref = d0_lin + alpha_fit .* v0_line;

    %% ── Forward model ODE ────────────────────────────────────────────────
    d_fwd   = nan(size(v0_line));
    r2_fwd  = NaN;

    if doForward
        dz        = 0.005;
        max_depth = max(d_all) * 2.0;
        v0_batch  = [v0_line(:); v0_all(:)];
        nB        = numel(v0_batch);
        v2        = v0_batch.^2;
        d_batch   = zeros(nB,1);
        stopped   = false(nB,1);
        z_now     = 0;

        while ~all(stopped) && z_now < max_depth
            F       = g_cm - k_over_m.*z_now - v2./d1;
            v2_new  = v2 + 2.*F.*dz;
            hit     = ~stopped & v2_new <= 0;
            d_batch(hit)  = z_now;
            stopped(hit)  = true;
            v2            = max(v2_new, 0);
            v2(stopped)   = 0;
            z_now         = z_now + dz;
        end
        d_batch(~stopped) = max_depth;

        nL            = numel(v0_line);
        d_fwd         = d_batch(1:nL)';
        d_pred_trials = d_batch(nL+1:end)';
        r2_fwd        = 1 - sum((d_all-d_pred_trials).^2)/ss_tot;

        fprintf('Forward model:  d1=%.3f cm  k/m=%.1f s^-2   R2=%.4f\n', ...
            d1, k_over_m, r2_fwd);
    end

    % 2. Depth-scaling with correct H(v0):
    %    h(v0) = v0^2/(2g), H = h + d
    %    Use d_fwd if available, else d_lin_ref as estimate of d
    if doForward
        d_est = d_fwd;
    else
        d_est = d_lin_ref;
    end
    h_line         = v0_line.^2 ./ (2*g_cm);
    H_line         = h_line + max(d_est, 0);
    d_scaling_line = (d0^2 .* H_line).^(1/3);

    fprintf('-------------------------------------------------------------\n\n');

    %% ── Figure ───────────────────────────────────────────────────────────
    fig = figure('Name','d vs v0','ToolBar','none','MenuBar','none');
    fig.Position = [100 100 680 540];
    ax = axes(fig, 'Position',[0.12 0.13 0.58 0.82]);
    hold(ax,'on');

    % 1. Linear reference — grey dashed
    plot(ax, v0_line, d_lin_ref, '--', 'Color',[0.35 0.35 0.35], ...
        'LineWidth',2.0, 'HandleVisibility','off');

    % 2. Depth-scaling — grey dash-dot (H corrected)
    plot(ax, v0_line, d_scaling_line, '-.', 'Color',[0.55 0.55 0.55], ...
        'LineWidth',2.0, 'HandleVisibility','off');

    % 3. Forward model — black dotted
    if doForward
        plot(ax, v0_line, d_fwd, 'k:', 'LineWidth',2.2, ...
            'HandleVisibility','off');
    end

    % Data points per height group
    for j = 1:numel(heights)
        hg = heights(j);
        [~, marker] = get_height_style(hg.h_cm);
        errorbar(ax, hg.v0_mean, hg.d_mean, hg.d_std, hg.d_std, ...
            hg.v0_std, hg.v0_std, marker, ...
            'Color',cmap(j,:), 'MarkerFaceColor','none', ...
            'MarkerEdgeColor',cmap(j,:), 'MarkerSize',9, 'LineWidth',1.8, ...
            'HandleVisibility','off');
    end

    set(ax,'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in', ...
        'XLim',[0, v0_max], ...
        'YLim',[0, max(d_all)*1.18]);
    grid(ax,'off');
    xlabel(ax,'$v_0$  (cm s$^{-1}$)','FontSize',16, ...
        'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax,'$d$  (cm)','FontSize',16, ...
        'Interpreter','latex','Color',[0 0 0]);

    %% ── Legend (dummy overlay axes) ──────────────────────────────────────
    ax2 = axes(fig,'Position',ax.Position,'Visible','off');
    hold(ax2,'on');

    plot(ax2, nan, nan, '--', 'Color',[0.35 0.35 0.35], 'LineWidth',2.0, ...
        'DisplayName', sprintf('$d = %.3f + %.5f\\,v_0$', d0_lin, alpha_fit));
    plot(ax2, nan, nan, '-.', 'Color',[0.55 0.55 0.55], 'LineWidth',2.0, ...
        'DisplayName', sprintf('$(d_0^2 H)^{1/3}$,  $d_0=%.3f$ cm', d0));
    if doForward
        plot(ax2, nan, nan, 'k:', 'LineWidth',2.2, ...
            'DisplayName', sprintf(...
                'Forward model ($d_1=%.2f$ cm, $k/m=%.0f$ s$^{-2}$)', ...
                d1, k_over_m));
    end

    lgd2 = legend(ax2,'show','FontSize',9,'Box','on', ...
        'Interpreter','latex','EdgeColor',[0.15 0.15 0.15],'Color',[1 1 1]);
    lgd2.Location = 'none';
    lgd2.Position = [0.72 0.60 0.26 0.17];

    %% ── Annotation ───────────────────────────────────────────────────────
    if doForward
        ann_str = sprintf(['$R^2_{\\rm lin} = %.3f$\n' ...
                           '$R^2_{\\rm scale} = %.3f$\n' ...
                           '$R^2_{\\rm fwd} = %.3f$\n' ...
                           '$d_0 = %.2f$ cm'], ...
            r2_lin, r2_scaling, r2_fwd, d0);
    else
        ann_str = sprintf(['$R^2_{\\rm lin} = %.3f$\n' ...
                           '$R^2_{\\rm scale} = %.3f$\n' ...
                           '$d_0 = %.2f$ cm'], ...
            r2_lin, r2_scaling, d0);
    end

    text(ax, 0.03, 0.97, ann_str, ...
        'FontSize',9,'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'Interpreter','latex','FontAngle','italic', ...
        'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25], ...
        'LineWidth',1.0,'Margin',4);
end