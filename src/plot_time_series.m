function [ax, peakVals] = plot_time_series(heights, field, scalarField, ...
                                            scaleFactor, figName, ylab, legendLoc, cmap)
% PLOT_TIME_SERIES  Raw dots + shaded mean ± 1 std band for any kinematics
%   field vs physical time.
%
%   Inputs:
%       heights     - struct array from group_trials_by_height
%       field       - kinematics field ('z_smooth','v_smooth','a_plus_g')
%       scalarField - scalar field for inset ('d_final_cm','v0_cm_s')
%       scaleFactor - multiply scalar by this (e.g. 10 for cm->mm)
%       figName     - figure window name string
%       ylab        - y-axis label string
%       legendLoc   - legend location string
%       cmap        - optional [nHeights x 3] colormap matrix (default: get_height_style)

    if nargin < 8, cmap = []; end

    fig = figure('Name', figName);
    fig.Position = [100 100 820 520];
    ax  = axes(fig); hold(ax, 'on');

    % ── Raw dots (faint black, per trial) ────────────────────────────────
    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials
            k = heights(j).trials(i).kinematics;
            plot(ax, k.t_s(k.impact_index+1:k.stopFrame), ...
                 k.(field)(k.impact_index+1:k.stopFrame), '.', ...
                 'Color',            [0 0 0 0.25], ...
                 'MarkerSize',       3, ...
                 'HandleVisibility', 'off');
        end
    end

    % ── Shaded mean ± 1 std band ─────────────────────────────────────────
    peakVals = struct();
    for j = 1:numel(heights)
        hg = heights(j);

        % Color: use cmap if provided, else get_height_style
        if ~isempty(cmap)
            col = cmap(j,:);
        else
            col = get_height_style(hg.h_cm);
        end

        [tPhys, yMean, yStd, nValid] = interp_to_norm_time(hg, field);

        % Thinner line if fewer than 3 trials
        lw = 2.5;
        if hg.nTrials < 3, lw = 1.2; end

        plot_shaded_band(ax, tPhys, yMean, yStd, nValid, col, ...
            sprintf('v_0 = %.0f \\pm %.0f cm s^{-1}', hg.v0_mean, hg.v0_std), ...
            1, lw);

        peakVals.(sprintf('h%d', round(hg.h_cm))) = struct( ...
            'h',    hg.h_cm, ...
            'v0',   hg.v0_mean, ...
            'v0std',hg.v0_std, ...
            'vals', hg.(scalarField) * scaleFactor);
    end

    % ── Axes formatting ───────────────────────────────────────────────────
    xline(ax, 0, '--k', 'HandleVisibility', 'off');
    legend(ax, 'show', 'Location', legendLoc, 'FontSize', 11, 'Box', 'off');
    set(ax, 'FontSize', 13, 'Box', 'off');
    grid(ax, 'on');
    xlabel(ax, 't  (s)', 'FontSize', 14);
    ylabel(ax, ylab,     'FontSize', 14);
end