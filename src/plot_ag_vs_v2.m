function [fig, d1_fit] = plot_ag_vs_v2(heights, cmap_depths, z_targets)
% PLOT_AG_VS_V2  Replicates Katsuragi & Durian 2007 Fig. 3a.
%   At each fixed depth zi, extract (v, a+g) from every trial.
%   Fit: a+g = F(zi)/m + v^2/d1 — shared d1, depth-specific intercepts.
%   Main plot: a+g vs v.  Inset: a+g vs v^2 (parallel lines = constant d1).
%
%   Inputs:
%       heights     - struct array from group_trials_by_height
%       cmap_depths - [numel(z_targets) x 3] colormap for depth levels
%       z_targets   - vector of fixed depths to sample (cm)
%
%   Outputs:
%       fig    - figure handle
%       d1_fit - fitted inertial drag length scale (cm)

    nZ     = numel(z_targets);
    tol    = 0.05;    % cm — depth tolerance window
    ag_cap = 12000;   % cm/s² — physical cap to remove noisy outliers

    % ── Collect (v, a+g) at each fixed depth ─────────────────────────────
    v_at_z  = cell(nZ, 1);
    ag_at_z = cell(nZ, 1);

    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials
            k   = heights(j).trials(i).kinematics;
            idx = k.impact_index : k.stopFrame;
            z   = k.z_smooth(idx);
            v   = k.v_smooth(idx);
            ag  = k.a_plus_g(idx);
            for zi = 1:nZ
                in = abs(z - z_targets(zi)) < tol & ...
                     isfinite(v) & isfinite(ag)   & ...
                     v > 0 & ag > 0 & ag < ag_cap;
                if any(in)
                    v_at_z{zi}  = [v_at_z{zi};  v(in)];
                    ag_at_z{zi} = [ag_at_z{zi}; ag(in)];
                end
            end
        end
    end

    % ── Global fit: shared d1, depth-specific F(zi)/m ────────────────────
    n_total = sum(cellfun(@numel, v_at_z));
    A_mat   = zeros(n_total, nZ + 1);
    b_vec   = zeros(n_total, 1);
    row     = 1;
    for zi = 1:nZ
        n_zi = numel(v_at_z{zi});
        if n_zi == 0, continue; end
        rows              = row : row + n_zi - 1;
        A_mat(rows, zi)   = 1;
        A_mat(rows, nZ+1) = v_at_z{zi}.^2;
        b_vec(rows)       = ag_at_z{zi};
        row               = row + n_zi;
    end
    coeffs    = lsqnonneg(A_mat, b_vec);
    Fz_over_m = coeffs(1:nZ);
    d1_fit    = 1 / coeffs(nZ+1);

    fprintf('\n── Fig 3a fit ───────────────────────────────────\n');
    fprintf('d1 = %.3f cm  (kinematic, from a+g vs v²)\n', d1_fit);
    for zi = 1:nZ
        fprintf('  z = %.2f cm:  F/m = %.1f cm s^{-2}\n', z_targets(zi), Fz_over_m(zi));
    end
    fprintf('─────────────────────────────────────────────────\n\n');

    % ── Figure ────────────────────────────────────────────────────────────
    fig = figure('Name','a+g vs v (Fig 3a)','ToolBar','none','MenuBar','none');
    fig.Position = [100 100 700 540];
    ax = axes(fig, 'Position', [0.12 0.13 0.55 0.82]);
    hold(ax, 'on');
    ag_max_plot = 0;

    for zi = 1:nZ
        if isempty(v_at_z{zi}), continue; end
        col      = cmap_depths(zi,:);
        v_max_zi = max(v_at_z{zi});
        v_ref    = linspace(0, v_max_zi*1.05, 200);
        ag_ref   = Fz_over_m(zi) + v_ref.^2 ./ d1_fit;
        ag_max_plot = max(ag_max_plot, max(ag_at_z{zi}));

        scatter(ax, v_at_z{zi}, ag_at_z{zi}, 30, col, 'filled', ...
            'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');
        plot(ax, v_ref, ag_ref, ':', 'Color', col, 'LineWidth', 2.2, ...
            'DisplayName', sprintf('z = %.2f cm', z_targets(zi)));
    end

    nonempty = ~cellfun(@isempty, v_at_z);
    set(ax, 'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in', ...
        'XLim',[0, max(cellfun(@max, v_at_z(nonempty)))*1.08], ...
        'YLim',[0, ag_max_plot*1.12]);
    grid(ax,'off');
    xlabel(ax,'$v$  (cm s$^{-1}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax,'$a+g$  (cm s$^{-2}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    lgd = legend(ax,'show','FontSize',11,'Box','on','Interpreter','tex', ...
        'EdgeColor',[0.25 0.25 0.25],'Color',[1 1 1]);
    lgd.Location = 'northwest';

    text(ax, 0.97, 0.05, sprintf('d_1 = %.2f cm', d1_fit), ...
        'FontSize',12,'Units','normalized', ...
        'HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'FontAngle','italic','Interpreter','tex','Color',[0.15 0.15 0.15], ...
        'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25], ...
        'LineWidth',1.0,'Margin',4);

    % ── Inset: a+g vs v² ──────────────────────────────────────────────────
    ax_in = axes(fig, 'Position', [0.72 0.58 0.24 0.34]);
    hold(ax_in, 'on');
    for zi = 1:nZ
        if isempty(v_at_z{zi}), continue; end
        col    = cmap_depths(zi,:);
        v2_max = max(v_at_z{zi}.^2);
        v2_ref = linspace(0, v2_max*1.05, 200);
        ag_ref = Fz_over_m(zi) + v2_ref ./ d1_fit;
        scatter(ax_in, v_at_z{zi}.^2, ag_at_z{zi}, 12, col, 'filled', ...
            'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none','HandleVisibility','off');
        plot(ax_in, v2_ref, ag_ref, ':', 'Color', col, 'LineWidth', 1.8, ...
            'HandleVisibility','off');
    end
    set(ax_in,'FontSize',10,'Box','on','LineWidth',1.0, ...
        'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','in');
    grid(ax_in,'off');
    xlabel(ax_in,'$v^2$  (cm$^2$ s$^{-2}$)','FontSize',11,'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax_in,'$a+g$','FontSize',11,'Interpreter','latex','Color',[0 0 0]);
end
