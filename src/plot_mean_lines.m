function ax = plot_mean_lines(heights, field, ylab, cmap)
% PLOT_MEAN_LINES  Mean lines only — no raw dots, no shaded band, no inset.
%   Uses interp_to_norm_time to compute mean across trials per height,
%   then rescales to physical time via mean t_stop.
%
%   Inputs:
%       heights - struct array from group_trials_by_height
%       field   - kinematics field ('z_smooth','v_smooth','a_plus_g')
%       ylab    - y-axis label string
%       cmap    - [nHeights x 3] colormap matrix
%
%   Output:
%       ax - axes handle for further customization

    fig = figure('Name', [field ' clean']);
    fig.Position = [100 100 820 520];
    ax  = axes(fig, 'Position', [0.10 0.12 0.72 0.82]);
    hold(ax, 'on');

    for j = 1:numel(heights)
        hg  = heights(j);
        col = cmap(j,:);
        [tPhys, yMean, ~, nValid] = interp_to_norm_time(hg, field);
        ok  = nValid >= 1 & isfinite(yMean);
        plot(ax, tPhys(ok), yMean(ok), '-', ...
            'Color',       col, ...
            'LineWidth',   2.5, ...
            'DisplayName', sprintf('v_0 = %.0f cm s^{-1}', hg.v0_mean));
    end
    xline(ax, 0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, ...
        'HandleVisibility', 'off');

    % Legend outside to the right
    lgd          = legend(ax, 'show', 'FontSize', 10, 'Box', 'off');
    lgd.Location = 'none';
    lgd.Position = [0.84 0.30 0.14 0.50];

    set(ax, 'FontSize', 13, 'Box', 'on', 'LineWidth', 1.2, ...
        'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);
    grid(ax, 'off');
    xlabel(ax, 't  (s)', 'FontSize', 14);
    ylabel(ax, ylab,     'FontSize', 14);
end
