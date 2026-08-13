function [fig_z, fig_v, fig_ag] = plot_example_trials(heights, cmap)
% PLOT_EXAMPLE_TRIALS  Show z, v, and a+g vs t for three representative
%   trials: lowest, intermediate, and highest mean v0 height group.
%   Panels are sized to be 1:1 (square) in pixel space.

    nH    = numel(heights);
    i_lo  = max(2, round(nH * 0.20));
    i_mid = round(nH/2);
    i_hi  = nH;
    sel   = [i_lo, i_mid, i_hi];

    fields   = {'z_smooth',  'v_smooth',  'a_plus_g'};
    ylabels  = {'$z$  (cm)', '$v$  (cm s$^{-1}$)', '$a+g$  (cm s$^{-2}$)'};
    figNames = {'z vs t — examples', 'v vs t — examples', 'a+g vs t — examples'};
    figs     = gobjects(3,1);

    % ── Layout: compute square panel dimensions ────────────────────────────
    fig_w    = 1080;
    l_margin = 0.09;
    r_margin = 0.02;
    b_margin = 0.14;
    t_margin = 0.10;
    gap      = 0.04;

    % Panel width (normalized) — same for all panels
    pw = (1 - l_margin - r_margin - 2*gap) / 3;   % ≈ 0.27

    % Figure height so that panel height = panel width in pixels
    pw_px = pw * fig_w;                            % panel width in pixels
    fig_h = ceil(pw_px / (1 - b_margin - t_margin));
    ph    = pw_px / fig_h;                         % panel height (normalized)

    % Panel left edges
    pl = l_margin + (0:2) .* (pw + gap);           % [p1_left, p2_left, p3_left]

    for fi = 1:3

        figs(fi) = figure('Name', figNames{fi}, ...
                          'ToolBar','none','MenuBar','none');
        figs(fi).Position = [80 + (fi-1)*20, 100, fig_w, fig_h];

        ax_handles = gobjects(3,1);

        for pi = 1:3
            hg  = heights(sel(pi));
            col = cmap(sel(pi),:);

            % Manual square axes
            ax = axes(figs(fi), 'Position', [pl(pi), b_margin, pw, ph]);
            ax_handles(pi) = ax;
            hold(ax,'on');

            % ── Raw dots ──────────────────────────────────────────────────
            for i = 1:hg.nTrials
                k        = hg.trials(i).kinematics;
                idx_plot = 1:5:k.stopFrame;
                plot(ax, k.t_s(idx_plot), k.(fields{fi})(idx_plot), '.', ...
                    'Color',            [col, 0.20], ...
                    'MarkerSize',       4, ...
                    'LineStyle',        'none', ...
                    'HandleVisibility', 'off');
            end

            % ── Shaded ±0.5σ band + mean line ─────────────────────────────
            [tPhys, yMean, yStd, nValid] = interp_to_norm_time(hg, fields{fi});
            lw = 2.5;
            if hg.nTrials < 3, lw = 1.2; end
            plot_shaded_band(ax, tPhys, yMean, yStd, nValid, col, '', 1, lw);

            % ── Reference lines ────────────────────────────────────────────
            if fi ~= 3
                xline(ax, 0, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.4, ...
                    'HandleVisibility', 'off');
            end

            if fi == 1
                yline(ax, 0, '--k', 'LineWidth', 1, ...
                    'HandleVisibility', 'off');

            elseif fi == 2
                yline(ax, 0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, ...
                    'HandleVisibility', 'off');
                [~, marker] = get_height_style(hg.h_cm);
                errorbar(ax, 0, hg.v0_mean, hg.v0_std, marker, ...
                    'Color', col, 'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', col, 'MarkerSize', 9, ...
                    'LineWidth', 1.8, 'HandleVisibility', 'off');
                ylim(ax, [0, ax.YLim(2)]);
            end

            % ── Axes ───────────────────────────────────────────────────────
            set(ax, 'FontSize',   13,  'Box',        'on', ...
                    'LineWidth',  1.2, 'XColor',     [0 0 0], ...
                    'YColor',     [0 0 0], ...
                    'XMinorTick', 'on', 'YMinorTick', 'on', ...
                    'TickDir',    'in');
            grid(ax, 'off');

            if fi == 3
                set(ax, 'XLim', [0, 0.045]);
                ylim(ax, [0, ax.YLim(2)]);
            else
                set(ax, 'XLim', [-0.005, 0.045]);
            end

            xlabel(ax, '$t$  (s)', 'FontSize', 16, ...
                'Interpreter', 'latex', 'Color', [0 0 0]);
            ylabel(ax, ylabels{fi}, 'FontSize', 16, ...
                'Interpreter', 'latex', 'Color', [0 0 0]);

            title(ax, sprintf('$v_0 = %.0f$ cm s$^{-1}$', hg.v0_mean), ...
                'Interpreter', 'latex', 'FontSize', 12, ...
                'FontWeight', 'bold', 'Color', [0 0 0]);
        end

        % Link x-axes across all three panels
        linkaxes(ax_handles, 'x');
    end

    fig_z  = figs(1);
    fig_v  = figs(2);
    fig_ag = figs(3);
end