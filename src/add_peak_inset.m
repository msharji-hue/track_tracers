function ax = add_peak_inset(peakVals, hVals, physMeasured, pos, xlab, ylab)
% ADD_PEAK_INSET  Inset with hollow markers — shape matched between
%   imaging (color, hollow) and physical (black, hollow, same shape).

    ax = axes('Position', pos);
    hold(ax, 'on');

    % ── Imaging markers (hollow, colored, per height) ─────────────────────
    fields = fieldnames(peakVals);
    for fi = 1:numel(fields)
        d      = peakVals.(fields{fi});
        vals   = d.vals(:);
        mu     = mean(vals, 'omitnan');
        sg     = std(vals,  0, 'omitnan');
        [col, marker] = get_height_style(d.h);

        errorbar(ax, d.h, mu, sg, marker, ...
            'Color',            col, ...
            'MarkerFaceColor',  'none', ...    % hollow
            'MarkerEdgeColor',  col, ...
            'MarkerSize',       9, ...
            'LineWidth',        1.8, ...
            'HandleVisibility', 'off');
    end

    % ── Physical measurements (hollow, black, same shape as height) ───────
    if ~isempty(physMeasured)
        for j = 1:size(physMeasured, 1)
            h_j = physMeasured(j, 1);
            mu_j = physMeasured(j, 2);
            sg_j = physMeasured(j, 3);
            if isnan(mu_j), continue; end

            [~, marker] = get_height_style(h_j);   % same shape as imaging

            errorbar(ax, h_j, mu_j, sg_j, marker, ...
                'Color',            [0 0 0], ...
                'MarkerFaceColor',  'none', ...    % hollow
                'MarkerEdgeColor',  [0 0 0], ...
                'MarkerSize',       9, ...
                'LineWidth',        1.5, ...
                'HandleVisibility', 'off');
        end
    end

    % ── Legend entries ────────────────────────────────────────────────────
    plot(ax, nan, nan, 'ks', ...
        'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'DisplayName', 'imaging');
    plot(ax, nan, nan, 'ko', ...
        'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'DisplayName', 'physical');
    legend(ax, 'show', 'FontSize', 8, 'Box', 'off', 'Location', 'northwest');

    % ── Formatting ────────────────────────────────────────────────────────
    set(ax, 'FontSize', 10, 'Box', 'on', 'LineWidth', 1);
    xlabel(ax, xlab, 'FontSize', 11);
    ylabel(ax, ylab, 'FontSize', 11);
    grid(ax, 'on');

    if ~isempty(hVals)
        rng = max(hVals) - min(hVals);
        pad = max(rng * 0.1, 2);
        set(ax, 'XLim', [min(hVals) - pad, max(hVals) + pad]);
    end
end