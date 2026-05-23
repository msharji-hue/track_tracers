function fig = plot_fz_vs_z(heights, d1, cmap)
% PLOT_FZ_VS_Z  Plot F(z)/m = a+g - v²/d1 vs z for all trials.
%
%   This is the direct test of the Katsuragi friction term k|z|.
%   If the force law F = -mg + k|z| + mv²/d1 is correct, then
%   F(z)/m = a+g - v²/d1 should collapse to a single line: k/m * z
%   regardless of drop height (i.e. regardless of v0).
%
%   A good collapse to a single line through the origin confirms:
%     (1) d1 is correctly estimated from the kinematic fit
%     (2) Coulomb friction is depth-linear (k|z|)
%     (3) Force law is separable: inertial + friction terms decouple
%
%   Deviation from linearity (e.g. super-linear growth) would suggest
%   a jamming volume V_j contribution as proposed by Moore et al.
%
%   Inputs:
%       heights - struct array from group_trials_by_height
%       d1      - inertial drag length scale (cm), from plot_ag_vs_v2
%       cmap    - [nH x 3] colormap

    nH     = numel(heights);
    ag_cap = 8000;   % cm/s² — discard noisy outliers near impact

    fig = figure('Name','F(z)/m vs z','ToolBar','none','MenuBar','none');
    fig.Position = [100 100 620 500];
    ax = axes(fig, 'Position', [0.14 0.13 0.76 0.82]);
    hold(ax, 'on');

    z_all  = [];
    fz_all = [];

    for j = 1:nH
        for i = 1:heights(j).nTrials
            k   = heights(j).trials(i).kinematics;
            idx = k.impact_index:k.stopFrame;
            z   = k.z_smooth(idx);
            v   = k.v_smooth(idx);
            ag  = k.a_plus_g(idx);

            fz  = ag - v.^2 ./ d1;   % F(z)/m with inertial term removed

            ok  = isfinite(z) & isfinite(fz) & z > 0 & ag < ag_cap & ag > 0;
            if ~any(ok), continue; end

            plot(ax, z(ok), fz(ok), '.', ...
                'Color', [cmap(j,:), 0.40], 'MarkerSize', 4, ...
                'HandleVisibility', 'off');

            z_all  = [z_all;  z(ok)];
            fz_all = [fz_all; fz(ok)];
        end
    end

    % ── Linear fit through origin: F(z)/m = (k/m)*z ──────────────────────
    if numel(z_all) > 10
        k_over_m = z_all(:) \ fz_all(:);   % OLS with no intercept
        z_ref    = linspace(0, max(z_all)*1.05, 200);
        fz_ref   = k_over_m .* z_ref;
        plot(ax, z_ref, fz_ref, 'k-', 'LineWidth', 2.0, ...
            'DisplayName', sprintf('$k/m = %.0f$ s$^{-2}$', k_over_m));

        % Annotation
        text(ax, 0.05, 0.92, ...
            sprintf('$k/m = %.0f$ cm$^{-1}$ s$^{-2}$\n$d_1 = %.2f$ cm', ...
                k_over_m, d1), ...
            'Units', 'normalized', 'Interpreter', 'latex', ...
            'FontSize', 11, 'VerticalAlignment', 'top', ...
            'Color', [0.10 0.10 0.10], ...
            'BackgroundColor', [1 1 1], 'EdgeColor', [0.25 0.25 0.25], ...
            'Margin', 4);

        fprintf('\n── F(z)/m vs z fit ──────────────────────────────\n');
        fprintf('k/m = %.1f cm^-1 s^-2  (used d1 = %.3f cm)\n', k_over_m, d1);
        fprintf('─────────────────────────────────────────────────\n\n');
    end

    set(ax, 'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in', ...
        'XLim',[0 max(z_all)*1.08], 'YLim',[0 max(fz_all(isfinite(fz_all)))*1.12]);
    grid(ax,'off');
    xlabel(ax,'$z$  (cm)',            'FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax,'$(a+g) - v^2/d_1$  (cm s$^{-2}$)', ...
        'FontSize',14,'Interpreter','latex','Color',[0 0 0]);

    legend(ax,'show','FontSize',11,'Interpreter','latex','Box','off','Location','northwest');

    % Collapse quality annotation
    text(ax, 0.97, 0.06, ...
        {'Collapse $\Rightarrow$ $k|z|$ confirmed', ...
         'Scatter $\Rightarrow$ jamming volume $V_j$ needed'}, ...
        'Units','normalized','Interpreter','latex','FontSize',10, ...
        'HorizontalAlignment','right','Color',[0.25 0.25 0.25]);
end
