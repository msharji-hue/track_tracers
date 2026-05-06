function fig = plot_normalized_collapse(heights, cmap)

    % ── Collect normalized data per trial ─────────────────────────────────
    z_all = []; v_all = []; trial_z = {}; trial_v = {}; ti = 0;

    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials
            k     = heights(j).trials(i).kinematics;
            v0    = k.v_smooth(k.impact_index);
            z_max = max(k.z_smooth, [], 'omitnan');
            if v0 == 0 || z_max == 0; continue; end
            z_n   = k.z_smooth ./ z_max;
            v_n   = k.v_smooth ./ v0;
            valid = isfinite(z_n) & isfinite(v_n) & z_n >= 0 & z_n <= 1;
            ti    = ti + 1;
            trial_z{ti} = z_n(valid);
            trial_v{ti} = v_n(valid);
            z_all = [z_all; z_n(valid)];
            v_all = [v_all; v_n(valid)];
        end
    end

    % ── Alpha fit ─────────────────────────────────────────────────────────
    mask      = z_all < 0.99 & v_all > 0;
    alpha_fit = log(1 - z_all(mask)) \ log(v_all(mask));

    % ── Trial-level bootstrap 90% CI ──────────────────────────────────────
    rng(42);
    alpha_boot = zeros(500, 1);
    for kb = 1:500
        idx_t          = randi(ti, ti, 1);
        z_s            = vertcat(trial_z{idx_t});
        v_s            = vertcat(trial_v{idx_t});
        m              = z_s < 0.99 & v_s > 0;
        alpha_boot(kb) = log(1 - z_s(m)) \ log(v_s(m));
    end
    alpha_lo = prctile(alpha_boot, 5);
    alpha_hi = prctile(alpha_boot, 95);

    % ── Interpolate each trial onto shared z-grid, then average ───────────
    zRef = linspace(0, 1, 300);
    nH   = numel(heights);
    v_mean = nan(nH, numel(zRef));
    v_std  = nan(nH, numel(zRef));

    for j = 1:nH
        vMat = [];
        for i = 1:heights(j).nTrials
            k     = heights(j).trials(i).kinematics;
            v0    = k.v_smooth(k.impact_index);
            z_max = max(k.z_smooth, [], 'omitnan');
            if v0 == 0 || z_max == 0; continue; end
            z_n = k.z_smooth ./ z_max;
            v_n = k.v_smooth ./ v0;
            valid = isfinite(z_n) & isfinite(v_n);
            if sum(valid) < 5; continue; end
            [z_u, ia] = unique(z_n(valid));
            v_u = v_n(valid); v_u = v_u(ia);
            vi  = interp1(z_u, v_u, zRef, 'pchip', NaN);
            vMat = [vMat; vi];
        end
        if isempty(vMat); continue; end
        v_mean(j,:) = mean(vMat, 1, 'omitnan');
        v_std(j,:)  = std(vMat,  0, 1, 'omitnan');
    end

    % ── Manually separate the two darkest blues ───────────────────────────
    cmap(1,:) = [0.10 0.10 0.55];
    cmap(2,:) = [0.35 0.55 0.95];

    % ── Figure ────────────────────────────────────────────────────────────
    fig = figure('Name', 'normalized collapse', ...
                 'ToolBar', 'none', 'MenuBar', 'none');
    fig.Position = [100 100 960 580];
    ax = axes(fig, 'Position', [0.11 0.13 0.65 0.82]);
    hold(ax, 'on');

    % ── Bootstrap CI band ─────────────────────────────────────────────────
    fill(ax, [zRef, fliplr(zRef)], ...
        [max(0,1-zRef).^alpha_hi, fliplr(max(0,1-zRef).^alpha_lo)], ...
        [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    % ── Per-height smooth mean lines ───────────────────────────────────────
    for j = 1:nH
        col  = cmap(j,:);
        good = isfinite(v_mean(j,:));
        plot(ax, zRef(good), v_mean(j,good), '-', ...
            'Color', col, 'LineWidth', 2.5, ...
            'DisplayName', sprintf('v_0 = %.0f cm s^{-1}', heights(j).v0_mean));
    end

    % ── Reference curves ──────────────────────────────────────────────────
        plot(ax, zRef, sqrt(max(0, 1-zRef)), 'k:', 'LineWidth', 1.8, ...
            'DisplayName', '(1-z/z_{max})^{0.50}  [const. drag]');
        plot(ax, zRef, max(0, 1-zRef).^(2/3), '-.', ...
            'Color', [0.40 0.40 0.40], 'LineWidth', 1.8, ...
            'DisplayName', '(1-z/z_{max})^{0.67}  [linear resist.]');
        plot(ax, zRef, max(0, 1-zRef).^1.0, ':', ...        
            'Color', [0.7 0.3 0.3], 'LineWidth', 1.6, ...
            'DisplayName', '(1-z/z_{max})^{1.00}  [quadratic drag]');
        plot(ax, zRef, max(0, 1-zRef).^alpha_fit, 'k--', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('(1-z/z_{max})^{%.2f}  [fit]', alpha_fit));

    % ── Legend ────────────────────────────────────────────────────────────
    lgd = legend(ax, 'show', 'FontSize', 10, 'Box', 'off', ...
                 'Interpreter', 'tex');
    lgd.Location = 'none';
    lgd.Position = [0.78 0.18 0.20 0.70];

    % ── Axes ──────────────────────────────────────────────────────────────
    set(ax, 'FontSize', 13, 'Box', 'on', 'LineWidth', 1.2, ...
        'XColor', [0 0 0], 'YColor', [0 0 0], ...
        'XLim', [0 1], 'YLim', [0 1.02], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickDir', 'in');
    grid(ax, 'off');
    xlabel(ax, '$z \,/\, z_\mathrm{max}$', ...
        'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);
    ylabel(ax, '$v \,/\, v_0$', ...
        'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);

    % ── Alpha annotation (lower-left, boxed) ──────────────────────────────
    text(ax, 0.05, 0.12, ...
        sprintf('$\\alpha = %.2f\\;[%.2f,\\,%.2f]$', alpha_fit, alpha_lo, alpha_hi), ...
        'FontSize', 12, 'Units', 'normalized', 'FontAngle', 'italic', ...
        'Interpreter', 'latex', 'Color', [0.15 0.15 0.15], ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'BackgroundColor', [1 1 1 0.75], 'EdgeColor', [0.7 0.7 0.7], ...
        'Margin', 3);

    exportgraphics(fig, 'fig_collapse.pdf', 'ContentType', 'vector');
end