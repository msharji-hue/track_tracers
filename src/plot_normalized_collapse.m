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

    % ── Global alpha fit ──────────────────────────────────────────────────
    mask      = z_all < 0.99 & v_all > 0;
    alpha_fit = log(1 - z_all(mask)) \ log(v_all(mask));

    % ── Trial-level bootstrap 90% CI (global) ─────────────────────────────
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

    % ── Per-height alpha fit + within-group bootstrap CI ──────────────────
    nH        = numel(heights);
    alpha_j   = nan(nH, 1);
    alpha_jlo = nan(nH, 1);
    alpha_jhi = nan(nH, 1);
    v0vals    = nan(nH, 1);

    for j = 1:nH
        v0vals(j) = heights(j).v0_mean;
        zj = []; vj = []; gj_z = {}; gj_v = {}; gc = 0;
        for i = 1:heights(j).nTrials
            k     = heights(j).trials(i).kinematics;
            v0    = k.v_smooth(k.impact_index);
            z_max = max(k.z_smooth, [], 'omitnan');
            if v0 == 0 || z_max == 0; continue; end
            z_n   = k.z_smooth ./ z_max;
            v_n   = k.v_smooth ./ v0;
            valid = isfinite(z_n) & isfinite(v_n) & z_n >= 0 & z_n < 0.99 & v_n > 0;
            if ~any(valid); continue; end
            zj = [zj; z_n(valid)];
            vj = [vj; v_n(valid)];
            gc = gc + 1;
            gj_z{gc} = z_n(valid);
            gj_v{gc} = v_n(valid);
        end
        if numel(zj) < 10; continue; end
        alpha_j(j) = log(1 - zj) \ log(vj);

        if gc >= 2
            ab = nan(300, 1);
            for kb = 1:300
                ib     = randi(gc, gc, 1);
                zb     = vertcat(gj_z{ib});
                vb     = vertcat(gj_v{ib});
                mb     = zb < 0.99 & vb > 0;
                if sum(mb) < 5; continue; end
                ab(kb) = log(1 - zb(mb)) \ log(vb(mb));
            end
            ab = ab(isfinite(ab));
            if numel(ab) > 10
                alpha_jlo(j) = prctile(ab, 5);
                alpha_jhi(j) = prctile(ab, 95);
            end
        end
    end

    % ── Log-linear fit: alpha = a + b*ln(v0) ──────────────────────────────
    ok          = isfinite(alpha_j) & isfinite(v0vals);
    X_ln        = [ones(sum(ok),1), log(v0vals(ok))];
    p_ln        = X_ln \ alpha_j(ok);
    v0ref       = linspace(min(v0vals(ok))*0.92, max(v0vals(ok))*1.08, 200);
    alpha_trend = p_ln(1) + p_ln(2) .* log(v0ref);

    % Diagnostics
    alpha_pred   = p_ln(1) + p_ln(2) .* log(v0vals(ok));
    ss_res       = sum((alpha_j(ok) - alpha_pred).^2);
    ss_tot       = sum((alpha_j(ok) - mean(alpha_j(ok))).^2);
    r2_loglin    = 1 - ss_res / ss_tot;
    v0_geomean   = exp(mean(log(v0vals(ok))));
    v0_cross_050 = exp((0.50 - p_ln(1)) / p_ln(2));
    v0_cross_067 = exp((0.67 - p_ln(1)) / p_ln(2));

    fprintf('\nLog-linear fit: alpha = %.3f + %.3f * ln(v0)  |  R2 = %.3f\n', ...
        p_ln(1), p_ln(2), r2_loglin);
    fprintf('Centered form:  alpha = %.3f + %.3f * ln(v0 / %.1f)\n', ...
        alpha_fit, p_ln(2), v0_geomean);
    fprintf('Crossover  alpha=0.50: v0 = %.1f cm/s\n', v0_cross_050);
    fprintf('Crossover  alpha=0.67: v0 = %.1f cm/s\n', v0_cross_067);

    % ── Interpolate onto shared z-grid for mean lines ─────────────────────
    zRef   = linspace(0, 1, 300);
    v_mean = nan(nH, numel(zRef));

    for j = 1:nH
        vMat = [];
        for i = 1:heights(j).nTrials
            k     = heights(j).trials(i).kinematics;
            v0    = k.v_smooth(k.impact_index);
            z_max = max(k.z_smooth, [], 'omitnan');
            if v0 == 0 || z_max == 0; continue; end
            z_n   = k.z_smooth ./ z_max;
            v_n   = k.v_smooth ./ v0;
            valid = isfinite(z_n) & isfinite(v_n);
            if sum(valid) < 5; continue; end
            [z_u, ia] = unique(z_n(valid));
            v_u = v_n(valid); v_u = v_u(ia);
            vMat = [vMat; interp1(z_u, v_u, zRef, 'pchip', NaN)];
        end
        if ~isempty(vMat)
            v_mean(j,:) = mean(vMat, 1, 'omitnan');
        end
    end

    % ── Separate the two darkest blues ────────────────────────────────────
    cmap(1,:) = [0.10 0.10 0.55];
    cmap(2,:) = [0.35 0.55 0.95];

    % ── Figure ────────────────────────────────────────────────────────────
    fig = figure('Name', 'normalized collapse', ...
                 'ToolBar', 'none', 'MenuBar', 'none');
    fig.Position = [100 100 960 580];
    ax = axes(fig, 'Position', [0.10 0.13 0.76 0.82]);
    hold(ax, 'on');

    % ── Bootstrap CI band ─────────────────────────────────────────────────
    fill(ax, [zRef, fliplr(zRef)], ...
        [max(0,1-zRef).^alpha_hi, fliplr(max(0,1-zRef).^alpha_lo)], ...
        [0.5 0.5 0.5], 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    % ── Per-height mean lines ─────────────────────────────────────────────
    for j = 1:nH
        good = isfinite(v_mean(j,:));
        plot(ax, zRef(good), v_mean(j,good), '-', ...
            'Color', cmap(j,:), 'LineWidth', 2.2, ...
            'HandleVisibility', 'off');
    end

    % ── Reference curves ──────────────────────────────────────────────────
    plot(ax, zRef, max(0,1-zRef).^1.0, ':', ...
        'Color', [0.05 0.45 0.45], 'LineWidth', 2.4, ...
        'DisplayName', '$(1-z/z_\mathrm{max})^{1.00}$  [quadratic drag]');
    plot(ax, zRef, sqrt(max(0,1-zRef)), 'k:', 'LineWidth', 2.8, ...
        'DisplayName', '$(1-z/z_\mathrm{max})^{0.50}$  [const. drag]');
    plot(ax, zRef, max(0,1-zRef).^(2/3), '-.', ...
        'Color', [0.35 0.35 0.35], 'LineWidth', 2.8, ...
        'DisplayName', '$(1-z/z_\mathrm{max})^{0.67}$  [linear resist.]');
    plot(ax, zRef, max(0,1-zRef).^alpha_fit, 'k--', 'LineWidth', 3.4, ...
        'DisplayName', sprintf('$(1-z/z_\\mathrm{max})^{%.2f}$  [fit]', alpha_fit));

    % Reference legend — top right
    lgd = legend(ax, 'show', 'FontSize', 10, 'Box', 'on', ...
                 'Interpreter', 'latex', ...
                 'EdgeColor', [0.25 0.25 0.25], 'Color', [1 1 1]);
    lgd.Location = 'none';
    lgd.Position = [0.57 0.60 0.29 0.22];

    % ── Main axes formatting ──────────────────────────────────────────────
    set(ax, 'FontSize', 13, 'Box', 'on', 'LineWidth', 1.2, ...
        'XColor', [0 0 0], 'YColor', [0 0 0], ...
        'XLim', [0 1], 'YLim', [0 1.02], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'in');
    grid(ax, 'off');
    xlabel(ax, '$z \,/\, z_\mathrm{max}$', ...
        'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);
    ylabel(ax, '$v \,/\, v_0$', ...
        'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);
