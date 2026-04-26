function add_peak_inset(peakVals, hVals, physMeasured, pos, xlab, ylab)

    ax = axes('Position', pos);
    hold(ax, 'on');

    hArr = []; dArr = [];
    for fn = fieldnames(peakVals)'
        d      = peakVals.(fn{1});
        mu     = mean(d.vals);
        sg     = std(d.vals);
        [col, marker] = get_height_style(d.h);
        errorbar(ax, d.h, mu, sg, marker, 'Color', col, 'MarkerFaceColor', col, ...
            'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
        hArr(end+1) = d.h;
        dArr(end+1) = mu;
    end

    % Physical measurements
    if ~isempty(physMeasured)
        for j = 1:size(physMeasured,1)
            errorbar(ax, physMeasured(j,1), physMeasured(j,2), physMeasured(j,3), ...
                'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
    end

    set(ax, 'FontSize', 11, 'Box', 'on', 'LineWidth', 1);
    xlabel(ax, xlab, 'FontSize', 12);
    ylabel(ax, ylab, 'FontSize', 12);
    grid(ax, 'on');
end