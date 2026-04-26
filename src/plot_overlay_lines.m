function peakVals = plot_overlay_lines(ax, repTrials, trials, hVals, field, scalerField, scaleFactor)
% PLOT_OVERLAY_LINES  Plot one representative line per height, collect peak values.
%
%   field       - kinematics field to plot on y axis (e.g. 'z_smooth')
%   scalerField - scalars field for peak value (e.g. 'd_final_cm')
%   scaleFactor - multiply scalar by this (e.g. 10 for cm->mm, 1 otherwise)

    if nargin < 7, scaleFactor = 1; end
    peakVals = struct();

    % Representative lines
    for i = 1:numel(repTrials)
        k   = repTrials(i).kinematics;
        h   = repTrials(i).h_cm;
        col = get_height_style(h);
        plot(ax, k.t_s, k.(field), 'Color', col, 'LineStyle', '-', ...
            'LineWidth', 2.5, 'HandleVisibility', 'off');
    end

    % Collect peak values from all trials
    for i = 1:numel(trials)
        hkey = sprintf('h%d', round(trials(i).h_cm));
        if ~isfield(peakVals, hkey)
            peakVals.(hkey).h    = trials(i).h_cm;
            peakVals.(hkey).vals = [];
        end
        peakVals.(hkey).vals(end+1) = trials(i).scalars.(scalerField) * scaleFactor;
    end

    % Legend entries
    for h = hVals
        col = get_height_style(h);
        plot(ax, nan, nan, 'Color', col, 'LineStyle', '-', 'LineWidth', 2.5, ...
            'DisplayName', sprintf('h = %.2f cm', h));
    end
end
