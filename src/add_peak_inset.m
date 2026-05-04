function ax = add_peak_inset(peakVals, pos, xlab, ylab)
% ADD_PEAK_INSET  Inset axes with hollow markers on v0 x-axis.
%   Bidirectional error bars: vertical = scalar std, horizontal = v0 std.
%   Marker shape matches corresponding height in main plot.
%
%   Inputs:
%       peakVals - struct from plot_time_series (h, v0, v0std, vals per height)
%       pos      - axes Position [x y w h] in normalized units
%       xlab     - x-axis label string
%       ylab     - y-axis label string

    ax = axes('Position', pos);
    hold(ax, 'on');

    fields = fieldnames(peakVals);
    for fi = 1:numel(fields)
        d    = peakVals.(fields{fi});
        vals = d.vals(:);
        mu   = mean(vals, 'omitnan');
        sg   = std(vals,  0, 'omitnan');
        [col, marker] = get_height_style(d.h);

        errorbar(ax, d.v0, mu, sg, sg, d.v0std, d.v0std, marker, ...
            'Color',           col, ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', col, ...
            'MarkerSize',      8, ...
            'LineWidth',       1.8, ...
            'HandleVisibility','off');
    end

    set(ax, 'FontSize', 10, 'Box', 'on', 'LineWidth', 1.2, ...
        'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);
    xlabel(ax, xlab, 'FontSize', 11);
    ylabel(ax, ylab, 'FontSize', 11);
    grid(ax, 'off');
end