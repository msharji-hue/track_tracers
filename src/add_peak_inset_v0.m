function ax = add_peak_inset_v0(peakVals, cmap, pos, xlab, ylab)
% ADD_PEAK_INSET_V0  Inset scatter with peak value vs mean v0.
%   X-axis = mean v0 per height with horizontal ±1 std error bars.
%   Y-axis = mean scalar value with vertical ±1 std error bars.
%   Markers are hollow, colored by height via cmap.
%
%   Inputs:
%       peakVals - struct from plot_time_series, one field per height
%                  each with: .h, .v0, .v0std, .vals
%       cmap     - [nHeights x 3] colormap matching heights order
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
        col  = cmap(fi,:);
        [~, marker] = get_height_style(d.h);

        % Bidirectional error bars: horizontal = v0 spread, vertical = value spread
        errorbar(ax, d.v0, mu, sg, sg, d.v0std, d.v0std, marker, ...
            'Color',           col,  ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', col,  ...
            'MarkerSize',      9,    ...
            'LineWidth',       1.8,  ...
            'HandleVisibility','off');
    end

    set(ax, 'FontSize', 10, 'Box', 'on', 'LineWidth', 1);
    xlabel(ax, xlab, 'FontSize', 11);
    ylabel(ax, ylab, 'FontSize', 11);
    grid(ax, 'on');
end
