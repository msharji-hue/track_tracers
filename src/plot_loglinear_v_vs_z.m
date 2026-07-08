function fig = plot_loglinear_v_vs_z(heights, cmap, d1, k_over_m)
% PLOT_LOGLINEAR_V_VS_Z  Log-linear plot of v(z) for all trials.
%
%   Tests whether velocity decay is exponential (pure inertial) or
%   deviates from exponential due to friction contribution.
%
%   On log-linear axes:
%     - Pure exponential decay: straight line, slope = -1/d1
%     - Curvature above the line: friction slowing the decay
%     - Steeper slope at shallow z: local d1 larger near surface
%
%   Shows:
%     - Measured v(z) trajectories colored by height group
%     - Pure inertial reference: v = v0*exp(-z/d1) per trial (grey dashed)
%     - Pure friction reference: v from v^2 = v0^2 + 2gz - (k/m)*z^2
%
%   Inputs:
%       heights   - struct array
%       cmap      - [nH x 3] colormap
%       d1        - shared inertial length scale [cm]
%       k_over_m  - friction coefficient [s^-2]

    g_cm = 980;

    fig = figure('Name','Log-linear v vs z','ToolBar','none','MenuBar','none');
    fig.Position = [100 100 760 580];
    ax = axes(fig,'Position',[0.12 0.12 0.83 0.83]);
    hold(ax,'on');
    set(ax,'YScale','log', ...
        'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in');

    % Collect z range for reference lines
    z_max_all = 0;
    v_min_all = Inf;
    v_max_all = 0;

    % Plot measured trajectories + pure inertial reference per trial
    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials
            k_t  = heights(j).trials(i).kinematics;
            s_t  = heights(j).trials(i).scalars;
            idx  = k_t.impact_index:k_t.stopFrame;
            z_m  = k_t.z_smooth(idx);
            v_m  = k_t.v_smooth(idx);
            v0_t = s_t.v0_cm_s;

            ok = isfinite(z_m) & isfinite(v_m) & v_m > 0 & z_m >= 0;
            if sum(ok) < 3, continue; end

            z_ok = z_m(ok);
            v_ok = v_m(ok);

            % Measured trajectory
            plot(ax, z_ok, v_ok, '-', ...
                'Color',[cmap(j,:) 0.25], 'LineWidth',0.9, ...
                'HandleVisibility','off');

            % Pure inertial reference for this trial
            % v = v0 * exp(-z/d1)
            z_ref = linspace(0, max(z_ok), 200);
            v_inertial = v0_t .* exp(-z_ref ./ d1);
            plot(ax, z_ref, v_inertial, '--', ...
                'Color',[0.65 0.65 0.65], 'LineWidth',0.8, ...
                'HandleVisibility','off');

            z_max_all = max(z_max_all, max(z_ok));
            v_min_all = min(v_min_all, min(v_ok));
            v_max_all = max(v_max_all, max(v_ok));
        end
    end

    % Pure friction reference — one curve per height group
    % v^2 = v0^2 + 2gz - (k/m)*z^2, stops when v^2 = 0
    z_line = linspace(0, z_max_all, 300);
    for j = 1:numel(heights)
        v0_j  = heights(j).v0_mean;
        v2_fr = v0_j^2 + 2*g_cm.*z_line - k_over_m.*z_line.^2;
        ok_fr = v2_fr > 0;
        v_fr  = sqrt(max(v2_fr, 0));
        if any(ok_fr)
            plot(ax, z_line(ok_fr), v_fr(ok_fr), ':', ...
                'Color',[cmap(j,:) 0.6], 'LineWidth',1.4, ...
                'HandleVisibility','off');
        end
    end

    % Dummy legend entries
    plot(ax, nan, nan, '-',  'Color',[0.40 0.40 0.40], 'LineWidth',1.5, ...
        'DisplayName','Measured v(z)');
    plot(ax, nan, nan, '--', 'Color',[0.65 0.65 0.65], 'LineWidth',1.2, ...
        'DisplayName',sprintf('Pure inertial:  v = v_0 exp(-z/d_1),  d_1 = %.2f cm', d1));
    plot(ax, nan, nan, ':',  'Color',[0.40 0.40 0.40], 'LineWidth',1.4, ...
        'DisplayName',sprintf('Pure friction:  v^2 = v_0^2 + 2gz - (k/m)z^2'));

    set(ax,'XLim',[0, z_max_all*1.08], ...
           'YLim',[max(1, v_min_all*0.8), v_max_all*1.15]);
    grid(ax,'off');
    xlabel(ax,'$z$  (cm)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax,'$v$  (cm s$^{-1}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    legend(ax,'show','FontSize',10,'Interpreter','none','Box','on', ...
        'EdgeColor',[0.25 0.25 0.25],'Location','northeast');

    % Annotation
    ann = sprintf(['Log-linear axes: straight line = pure exponential\n' ...
                   'Curvature above line = friction contribution\n' ...
                   'd_1 = %.3f cm  (slope of inertial reference)\n' ...
                   'k/m = %.0f s^{-2}'], d1, k_over_m);
    text(ax, 0.97, 0.05, ann, ...
        'Units','normalized','Interpreter','none','FontSize',11, ...
        'Color',[0.08 0.08 0.08], ...
        'HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25],'Margin',5);
end
