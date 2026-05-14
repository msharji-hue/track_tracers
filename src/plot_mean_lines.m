function ax = plot_mean_lines(heights, field, ylab, cmap, varargin)
% PLOT_MEAN_LINES  Mean lines only — no raw dots, no shaded band, no inset.
%   Optional name-value args:
%       'xline'   true/false  — draw dashed vertical at t=0 (default true)
%       'yline'   true/false  — draw dashed horizontal at y=0 (default false)
%       'v0marks' true/false  — draw v0 errorbar markers at t=0 (default false)

    p = inputParser;
    addParameter(p, 'xline',   true,  @islogical);
    addParameter(p, 'yline',   false, @islogical);
    addParameter(p, 'v0marks', false, @islogical);
    addParameter(p, 'xlim0',   false, @islogical);  % clip x at 0
    addParameter(p, 'ylim0',   false, @islogical);  % clip y at 0
    parse(p, varargin{:});

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

    if p.Results.xline
        xline(ax, 0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end
    if p.Results.yline
        yline(ax, 0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end
    if p.Results.v0marks
        for j = 1:numel(heights)
            hg          = heights(j);
            [~, marker] = get_height_style(hg.h_cm);
            errorbar(ax, 0, hg.v0_mean, hg.v0_std, marker, ...
                'Color',           cmap(j,:), ...
                'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', cmap(j,:), ...
                'MarkerSize',      9, ...
                'LineWidth',       1.8, ...
                'HandleVisibility','off');
        end
    end

    if p.Results.xlim0, set(ax, 'XLim', [0 ax.XLim(2)]); end
    if p.Results.ylim0, set(ax, 'YLim', [0 ax.YLim(2)]); end

    lgd          = legend(ax, 'show', 'FontSize', 10, 'Box', 'off');
    lgd.Location = 'none';
    lgd.Position = [0.84 0.30 0.14 0.50];

    set(ax, 'FontSize',   13,  'Box',        'on', ...
            'LineWidth',  1.2, 'XColor',     [0 0 0], ...
            'YColor',     [0 0 0], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'TickDir',    'in');
    grid(ax, 'off');
    xlabel(ax, '$t$  (s)', 'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);
    ylabel(ax, ylab,       'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);
end