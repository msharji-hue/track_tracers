function alpha_fit = plot_normalized_collapse(ax, repTrials, nOverlay, hVals)
% PLOT_NORMALIZED_COLLAPSE  Plot v/v0 vs z/zmax with bootstrap CI and power law fit.

    % Collect normalized data
    z_all = []; v_all = []; trial_z = cell(nOverlay,1); trial_v = cell(nOverlay,1);
    for i = 1:nOverlay
        k     = repTrials(i).kinematics;
        v0    = k.v_smooth(k.impact_index);
        z_max = max(k.z_smooth, [], 'omitnan');
        if v0 == 0 || z_max == 0, continue; end
        z_n = k.z_smooth ./ z_max;  v_n = k.v_smooth ./ v0;
        valid = isfinite(z_n) & isfinite(v_n) & z_n >= 0 & z_n <= 1;
        trial_z{i} = z_n(valid);  trial_v{i} = v_n(valid);
        z_all = [z_all; z_n(valid)];  v_all = [v_all; v_n(valid)];
    end

    mask      = z_all < 0.99 & v_all > 0;
    alpha_fit = log(1 - z_all(mask)) \ log(v_all(mask));

    % Bootstrap 90% CI
    rng(42); n_boot = 500; alpha_boot = zeros(n_boot,1);
    for kb = 1:n_boot
        idx_t = randi(nOverlay, nOverlay, 1);
        z_s   = vertcat(trial_z{idx_t});  v_s = vertcat(trial_v{idx_t});
        m     = z_s < 0.99 & v_s > 0;
        alpha_boot(kb) = log(1 - z_s(m)) \ log(v_s(m));
    end
    alpha_lo = prctile(alpha_boot, 5);
    alpha_hi = prctile(alpha_boot, 95);

    zRef = linspace(0, 1, 300);
    fill(ax, [zRef, fliplr(zRef)], [max(0,1-zRef).^alpha_hi, fliplr(max(0,1-zRef).^alpha_lo)], ...
        [0.3 0.3 0.3], 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Data lines
    for i = 1:nOverlay
        k     = repTrials(i).kinematics;
        v0    = k.v_smooth(k.impact_index);
        z_max = max(k.z_smooth, [], 'omitnan');
        if v0 == 0 || z_max == 0, continue; end
        col   = get_height_style(repTrials(i).h_cm);
        plot(ax, k.z_smooth./z_max, k.v_smooth./v0, '-', 'Color', col, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end

    % Reference curves
    plot(ax, zRef, sqrt(max(0,1-zRef)), 'k:', 'LineWidth', 1.5, 'DisplayName', '\surd(1-z/z_{max})  [const. drag]');
    plot(ax, zRef, max(0,1-zRef).^alpha_fit, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 2.0, ...
        'DisplayName', sprintf('(1-z/z_{max})^{%.2f}  [fit]', alpha_fit));

    for h = hVals
        col = get_height_style(h);
        plot(ax, nan, nan, '-', 'Color', col, 'LineWidth', 2.5, 'DisplayName', sprintf('h = %.2f cm', h));
    end

    legend(ax, 'show', 'FontSize', 12, 'Box', 'off', 'Location', 'northeast');
    set(ax, 'FontSize', 13, 'Box', 'off', 'XLim', [0 1], 'YLim', [0 1.05]);
    xlabel(ax, 'z / z_{max}', 'FontSize', 14); ylabel(ax, 'v / v_0', 'FontSize', 14);
    grid(ax, 'on');
    text(ax, 0.97, 0.50, sprintf('\\alpha = %.2f  [%.2f, %.2f]', alpha_fit, alpha_lo, alpha_hi), ...
        'FontSize', 12, 'Units', 'normalized', 'FontAngle', 'italic', ...
        'Color', [0.2 0.2 0.2], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
end
