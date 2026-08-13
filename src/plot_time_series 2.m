function [ax, peakVals] = plot_time_series(heights, field, scalarField, ...
                                            scaleFactor, figName, ylab, legendLoc, cmap)
% PLOT_TIME_SERIES  Raw dots + shaded mean ± std band for any kinematics
%   field vs physical time. Used for z vs t, v vs t, a+g vs t overlays.
%
%   Inputs:
%       heights     - struct array from group_trials_by_height
%       field       - kinematics field ('z_smooth','v_smooth','a_plus_g')
%       scalarField - scalar field for inset ('d_final_cm','v0_cm_s')
%       scaleFactor - multiply scalar by this (e.g. 10 for cm->mm)
%       figName     - figure window name string
%       ylab        - y-axis label string
%       legendLoc   - 'none' to suppress, or any valid Location string
%       cmap        - optional [nHeights x 3] colormap (default: get_height_style)
%
%   Outputs:
%       ax       - axes handle
%       peakVals - struct with scalar values + v0 stats per height

    if nargin < 8, cmap = []; end

    fig = figure('Name', figName);
    fig.Position = [100 100 980 520];
    ax  = axes(fig, 'Position', [0.08 0.11 0.74 0.82]);
    hold(ax, 'on');

    % ── Raw dots: uniform '.' colored by height group, low opacity ────────
    for j = 1:numel(heights)

        if ~isempty(cmap)
            col = cmap(j,:);
        else
            col = get_height_style(heights(j).h_cm);
        end

        for i = 1:heights(j).nTrials
            k        = heights(j).trials(i).kinematics;
            idx_plot = 1:5:k.stopFrame;
            plot(ax, k.t_s(idx_plot), k.(field)(idx_plot), '.', ...
                'Color',            [col, 0.20], ...
                'MarkerSize',       4, ...
                'LineStyle',        'none', ...
                'HandleVisibility', 'off');
        end
    end

    % ── Shaded mean ± 0.5 std band (symmetric) ────────────────────────────
    peakVals = struct();

    for j = 1:numel(heights)
        hg = heights(j);

        if ~isempty(cmap)
            col = cmap(j,:);
        else
            col = get_height_style(hg.h_cm);
        end

        [tPhys, yMean, yStd, nValid] = interp_to_norm_time(hg, field);

        lw = 2.5;
        if hg.nTrials < 3, lw = 1.2; end

        plot_shaded_band(ax, tPhys, yMean, yStd, nValid, col, ...
            sprintf('v_0 = %.0f \\pm %.0f cm s^{-1}', hg.v0_mean, hg.v0_std), ...
            1, lw);

        peakVals.(sprintf('h%d', round(hg.h_cm))) = struct( ...
            'h',     hg.h_cm, ...
            'v0',    hg.v0_mean, ...
            'v0std', hg.v0_std, ...
            'vals',  hg.(scalarField) * scaleFactor);
    end

    % ── Axes formatting ───────────────────────────────────────────────────
    xline(ax, 0, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.4, ...
        'HandleVisibility', 'off');

    if ~strcmp(legendLoc, 'none')
        legend(ax, 'show', 'Location', legendLoc, 'FontSize', 11, 'Box', 'off');
    else
        lgd          = legend(ax, 'show', 'FontSize', 10, 'Box', 'off');
        lgd.Location = 'none';
        lgd.Position = [0.84 0.28 0.14 0.50];
    end

    set(ax, 'FontSize',   13,  'Box',        'on', ...
            'LineWidth',  1.2, 'XColor',     [0 0 0], ...
            'YColor',     [0 0 0], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'TickDir',    'in');
    grid(ax, 'off');
    xlabel(ax, '$t$  (s)', 'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);
    ylabel(ax, ylab,       'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);
end