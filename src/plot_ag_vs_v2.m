function [fig, d1_fit, Fz_over_m, fitStats] = plot_ag_vs_v2(heights, ~, ...
        z_targets_all, z_display, v_min)
% PLOT_AG_VS_V2  Katsuragi & Durian (2007) Fig. 3a — exact approach.
%
%   Two-stage regression:
%     Stage 1 — STABLE depths only (z >= z_min_fit):
%               Independent per-depth OLS of a+g vs v² at each depth.
%               Shared d1 = 1/mean(slopes) — Katsuragi's parallel-lines test.
%     Stage 2 — ALL 22 depths:
%               Per-depth OLS with d1 fixed from Stage 1.
%               Intercepts at shallow depths use local slope (not shared d1)
%               to avoid the negative-intercept artifact.
%
%   Returns Fz_over_m for all 22 depths — used by plot_fz_vs_z,
%   plot_fz_vs_V, and plot_fz_vs_Az.
%
%   Inputs:
%       heights       - struct array from group_trials_by_height
%       ~             - colormap (unused)
%       z_targets_all - all computation depths [cm], e.g. 0.10:0.10:2.2
%       z_display     - 5 depths to show in figure [cm]
%       v_min         - minimum velocity [cm/s], default 40
%
%   Outputs:
%       fig       - figure handle
%       d1_fit    - shared d1 from stable-depth slopes [cm]
%       Fz_over_m - [nZ_all x 1] OLS intercepts for ALL depths [cm s^-2]
%       fitStats  - struct with full per-depth details

    %% ── Defaults ─────────────────────────────────────────────────────────
    if nargin < 5, v_min = 40; end

    z_targets_all = z_targets_all(:);
    z_display     = z_display(:);
    nZ_all        = numel(z_targets_all);
    z_min_fit     = 1.30;   % stable shaft boundary
    tol           = 0.05;
    ag_cap        = 60000;
    nBins         = 16;

    % Display depth indices
    disp_idx = zeros(numel(z_display), 1);
    for di = 1:numel(z_display)
        [~, disp_idx(di)] = min(abs(z_targets_all - z_display(di)));
    end
    nZ_disp = numel(disp_idx);

    % Colors for displayed depths
    depth_cols = [0.05 0.05 0.05;
                  0.55 0.10 0.25;
                  0.88 0.25 0.55;
                  0.98 0.45 0.20;
                  0.20 0.45 0.80];
    if nZ_disp ~= 5
        depth_cols = interp1(linspace(0,1,5), depth_cols, linspace(0,1,nZ_disp));
    end

    %% ── STEP 1: Collect (v, a+g) at ALL depth targets ───────────────────
    v_raw  = cell(nZ_all, 1);
    ag_raw = cell(nZ_all, 1);

    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials
            k   = heights(j).trials(i).kinematics;
            idx = k.impact_index:k.stopFrame;
            z   = k.z_smooth(idx);
            v   = k.v_smooth(idx);
            ag  = k.a_plus_g(idx);

            for zi = 1:nZ_all
                in = abs(z - z_targets_all(zi)) < tol & ...
                     isfinite(v) & isfinite(ag) & ...
                     v > 0 & ag > 0 & ag < ag_cap;
                if any(in)
                    v_raw{zi}  = [v_raw{zi};  v(in)];
                    ag_raw{zi} = [ag_raw{zi}; ag(in)];
                end
            end
        end
    end

    % Apply v_min
    v_at_z  = cell(nZ_all, 1);
    ag_at_z = cell(nZ_all, 1);
    for zi = 1:nZ_all
        ok          = v_raw{zi} >= v_min;
        v_at_z{zi}  = v_raw{zi}(ok);
        ag_at_z{zi} = ag_raw{zi}(ok);
    end

    fprintf('\nplot_ag_vs_v2: data collection (ag_cap=%d, v_min=%.0f)\n', ...
        ag_cap, v_min);
    for zi = 1:nZ_all
        fprintf('  z=%.2f cm: %d raw -> %d after v_min\n', ...
            z_targets_all(zi), numel(v_raw{zi}), numel(v_at_z{zi}));
    end

    %% ── STEP 2: Per-depth OLS at ALL depths ──────────────────────────────
    % Each depth gets independent fit: a+g = b0 + b1*v²
    %   b0 = F(zi)/m  (intercept — Fig 3b symbol)
    %   b1 = 1/d1_local (slope)
    %
    % Stable depths (z >= z_min_fit): used to calibrate shared d1
    % All depths: intercepts returned for Step 3 figures

    b0_all    = nan(nZ_all, 1);   % intercepts F(zi)/m
    b1_all    = nan(nZ_all, 1);   % slopes 1/d1_local
    d1_local  = nan(nZ_all, 1);   % local d1 per depth
    r2_local  = nan(nZ_all, 1);   % R² per depth
    n_local   = zeros(nZ_all, 1); % n points per depth

    fprintf('\n-- Per-depth OLS at ALL depths (ag_cap=%d) ------------------\n', ag_cap);
    fprintf('  %-8s  %-6s  %-10s  %-10s  %-6s  %-6s\n', ...
        'z (cm)', 'n', 'F/m', 'd1_local', 'R²', 'regime');

    for zi = 1:nZ_all
        v2  = v_at_z{zi}.^2;
        ag  = ag_at_z{zi};
        n_local(zi) = numel(ag);

        if n_local(zi) < 3
            fprintf('  z=%.2f    %3d  [insufficient]\n', ...
                z_targets_all(zi), n_local(zi));
            continue
        end

        A_ols      = [ones(n_local(zi),1), v2];
        coeffs     = A_ols \ ag;
        b0_all(zi) = coeffs(1);
        b1_all(zi) = coeffs(2);
        d1_local(zi) = 1 / max(coeffs(2), 1e-10);

        ag_pred     = A_ols * coeffs;
        ss_res      = sum((ag - ag_pred).^2);
        ss_tot      = sum((ag - mean(ag)).^2);
        r2_local(zi)= 1 - ss_res/ss_tot;

        if z_targets_all(zi) >= z_min_fit
            regime = 'stable';
        elseif b0_all(zi) >= 0
            regime = 'transition+';
        else
            regime = 'transition-';
        end

        fprintf('  z=%.2f    %3d  %10.1f  %10.3f  %.3f  %s\n', ...
            z_targets_all(zi), n_local(zi), b0_all(zi), ...
            d1_local(zi), r2_local(zi), regime);
    end

    %% ── STEP 3: Shared d1 from stable depths only ────────────────────────
    stable_mask = z_targets_all >= z_min_fit & isfinite(b1_all) & b1_all > 0;
    ok_slopes   = stable_mask;
    d1_fit      = 1 / mean(b1_all(ok_slopes));
    d1_std      = std(d1_local(ok_slopes));
    d1_cv       = d1_std / d1_fit * 100;

    fprintf('\nShared d1 = %.3f cm  (mean of %d stable slopes)\n', ...
        d1_fit, sum(ok_slopes));
    fprintf('d1 std    = %.3f cm  (CV = %.1f%%)\n', d1_std, d1_cv);
    if d1_cv < 15
        fprintf('Slopes approximately parallel (CV < 15%%) — shared d1 justified\n');
    else
        fprintf('Notable slope variation (CV >= 15%%) — shared d1 approximate\n');
    end

    %% ── STEP 4: Fz_over_m — use OLS intercepts directly ─────────────────
    % For all depths: use the per-depth OLS intercept b0
    % This avoids the negative-intercept artifact of subtracting shared d1
    % at depths where local d1 differs significantly from shared d1
    Fz_over_m = b0_all;   % [nZ_all x 1], may be negative at very shallow z

    fprintf('\n-- Intercepts (all depths) -----------------------------------\n');
    for zi = 1:nZ_all
        if isfinite(Fz_over_m(zi))
            fprintf('  z=%.2f  F/m=%8.1f  d1_local=%.3f\n', ...
                z_targets_all(zi), Fz_over_m(zi), d1_local(zi));
        end
    end

    %% ── STEP 5: Per-height intercepts for stable depths ──────────────────
    nH               = numel(heights);
    Fz_per_height    = nan(nH, nZ_all);
    Fz_per_height_se = nan(nH, nZ_all);
    v0_per_height    = nan(nH, 1);

    for j = 1:nH
        v0_per_height(j) = heights(j).v0_mean;

        for zi = 1:nZ_all
            if ~stable_mask(zi), continue; end   % only stable depths
            v_hj = []; ag_hj = [];

            for i = 1:heights(j).nTrials
                k_h  = heights(j).trials(i).kinematics;
                idx  = k_h.impact_index:k_h.stopFrame;
                z_h  = k_h.z_smooth(idx);
                v_h  = k_h.v_smooth(idx);
                ag_h = k_h.a_plus_g(idx);

                in = abs(z_h - z_targets_all(zi)) < tol & ...
                     v_h >= v_min & isfinite(v_h) & isfinite(ag_h) & ...
                     v_h > 0 & ag_h > 0 & ag_h < ag_cap;

                v_hj  = [v_hj;  v_h(in)];
                ag_hj = [ag_hj; ag_h(in)];
            end

            if numel(v_hj) >= 3
                fz_vals               = ag_hj - v_hj.^2 ./ d1_fit;
                Fz_per_height(j,zi)   = mean(fz_vals);
                Fz_per_height_se(j,zi)= std(fz_vals) / sqrt(numel(fz_vals));
            end
        end
    end

    %% ── STEP 6: Pack fitStats ────────────────────────────────────────────
    fitStats.d1_fit           = d1_fit;
    fitStats.d1_local         = d1_local;
    fitStats.d1_std           = d1_std;
    fitStats.d1_cv            = d1_cv;
    fitStats.b0_all           = b0_all;        % all 22 OLS intercepts
    fitStats.b1_all           = b1_all;        % all 22 OLS slopes
    fitStats.r2_local         = r2_local;
    fitStats.n_local          = n_local;
    fitStats.stable_mask      = stable_mask;
    fitStats.Fz_over_m        = Fz_over_m;     % all 22, may include negatives
    fitStats.Fz_per_height    = Fz_per_height;
    fitStats.Fz_per_height_se = Fz_per_height_se;
    fitStats.v0_per_height    = v0_per_height;
    fitStats.v_min            = v_min;
    fitStats.z_min_fit        = z_min_fit;
    fitStats.z_targets_all    = z_targets_all;
    fitStats.z_display        = z_display;
    fitStats.disp_idx         = disp_idx;

    %% ── STEP 7: Figure — display subset only ────────────────────────────
    nonempty_disp = ~cellfun(@isempty, v_at_z(disp_idx));
    if ~any(nonempty_disp)
        error('plot_ag_vs_v2: no data at any display depth after v_min cutoff.');
    end

    v_disp_max = max(cellfun(@(c) max(c), v_at_z(disp_idx(nonempty_disp))));
    v_edges    = linspace(v_min, v_disp_max * 1.02, nBins + 1);
    v_centers  = (v_edges(1:end-1) + v_edges(2:end)) / 2;

    bin_mean = nan(nZ_disp, nBins);
    bin_std  = nan(nZ_disp, nBins);

    for di = 1:nZ_disp
        zi = disp_idx(di);
        if isempty(v_at_z{zi}), continue; end
        for b = 1:nBins
            in_b = v_at_z{zi} >= v_edges(b) & v_at_z{zi} < v_edges(b+1);
            if sum(in_b) >= 3
                bin_mean(di,b) = mean(ag_at_z{zi}(in_b));
                bin_std(di,b)  = std( ag_at_z{zi}(in_b));
            end
        end
    end

    ag_max_plot = max(bin_mean(:), [], 'omitnan') * 1.20;

    fig = figure('Name','a+g vs v (Katsuragi Fig 3a)', ...
                 'ToolBar','none','MenuBar','none');
    fig.Position = [100 100 740 560];
    ax = axes(fig, 'Position', [0.11 0.12 0.82 0.83]);
    hold(ax, 'on');

    v_line = linspace(0, v_disp_max * 1.04, 400);

    for di = 1:nZ_disp
        zi  = disp_idx(di);
        col = depth_cols(di,:);
        if isempty(v_at_z{zi}), continue; end

        % Faint raw scatter
        scatter(ax, v_at_z{zi}, ag_at_z{zi}, 5, col, 'filled', ...
            'MarkerFaceAlpha', 0.10, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % Binned means ± std
        ok = isfinite(bin_mean(di,:));
        errorbar(ax, v_centers(ok), bin_mean(di,ok), bin_std(di,ok), ...
            'o', 'Color', col, 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'w', 'MarkerSize', 9, ...
            'LineWidth', 1.4, 'CapSize', 3, ...
            'DisplayName', sprintf('$z = %.2f$ cm', z_targets_all(zi)));

        % Per-depth OLS fit line
        if isfinite(b0_all(zi)) && isfinite(b1_all(zi))
            ag_fit = b0_all(zi) + b1_all(zi) .* v_line.^2;
            plot(ax, v_line, ag_fit, '--', 'Color', col, 'LineWidth', 1.8, ...
                'HandleVisibility', 'off');
        end
    end

    tick_start = ceil(v_min / 50) * 50;
    set(ax, 'FontSize',14, 'Box','on', 'LineWidth',1.2, ...
        'XColor',[0 0 0], 'YColor',[0 0 0], ...
        'XMinorTick','on', 'YMinorTick','on', 'TickDir','in', ...
        'XLim',  [v_min * 0.9, v_disp_max * 1.05], ...
        'XTick', tick_start : 50 : ceil(v_disp_max/50)*50, ...
        'YLim',  [0, ag_max_plot]);
    grid(ax, 'off');
    xlabel(ax, '$v$  (cm s$^{-1}$)', ...
        'FontSize',20, 'Interpreter','latex', 'Color',[0 0 0]);
    ylabel(ax, '$a+g$  (cm s$^{-2}$)', ...
        'FontSize',20, 'Interpreter','latex', 'Color',[0 0 0]);
    legend(ax, 'show', 'FontSize',11, 'Box','on', 'Interpreter','latex', ...
        'EdgeColor',[0.25 0.25 0.25], 'Color',[1 1 1], 'Location','southeast');

    text(ax, 0.97, 0.05, ...
        sprintf(['$d_1 = %.2f$ cm\n' ...
                 'CV $= %.1f\\%%$ (%d stable depths)'], ...
            d1_fit, d1_cv, sum(ok_slopes)), ...
        'FontSize',10, 'Units','normalized', ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
        'FontAngle','italic', 'Interpreter','latex', ...
        'Color',[0.15 0.15 0.15], ...
        'BackgroundColor',[1 1 1], 'EdgeColor',[0.25 0.25 0.25], ...
        'LineWidth',1.0, 'Margin',4);

    % ── Inset: a+g vs v² ─────────────────────────────────────────────────
    ax_in = axes(fig, 'Position', [0.16 0.63 0.25 0.28]);
    hold(ax_in, 'on');

    v2_max     = v_disp_max^2 * 1.02;
    v2_edges   = linspace(v_min^2, v2_max, nBins+1);
    v2_centers = (v2_edges(1:end-1) + v2_edges(2:end)) / 2;

    for di = 1:nZ_disp
        zi  = disp_idx(di);
        col = depth_cols(di,:);
        if isempty(v_at_z{zi}), continue; end

        v2  = v_at_z{zi}.^2;
        bm2 = nan(1, nBins);
        for b = 1:nBins
            in_b = v2 >= v2_edges(b) & v2 < v2_edges(b+1);
            if sum(in_b) >= 3
                bm2(b) = mean(ag_at_z{zi}(in_b));
            end
        end
        ok2 = isfinite(bm2);
        plot(ax_in, v2_centers(ok2), bm2(ok2), 'o', ...
            'Color', col, 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 6, ...
            'HandleVisibility', 'off');

        if isfinite(b0_all(zi)) && isfinite(b1_all(zi))
            v2_line = linspace(0, v2_max, 300);
            plot(ax_in, v2_line, b0_all(zi) + b1_all(zi).*v2_line, ...
                '--', 'Color', col, 'LineWidth', 1.4, ...
                'HandleVisibility', 'off');
        end
    end

    set(ax_in, 'FontSize',10, 'Box','on', 'LineWidth',1.0, ...
        'XColor',[0 0 0], 'YColor',[0 0 0], 'TickDir','in', ...
        'XMinorTick','on', 'YMinorTick','on', ...
        'XLim',[0, v2_max], ...
        'YLim',[0, max(bin_mean(:),[],'omitnan')*1.15]);
    grid(ax_in, 'off');
    xlabel(ax_in, '$v^2$  (cm$^2$ s$^{-2}$)', ...
        'FontSize',15, 'Interpreter','latex', 'Color',[0 0 0]);
    ylabel(ax_in, '$a+g$  (cm s$^{-2}$)', ...
        'FontSize',15, 'Interpreter','latex', 'Color',[0 0 0]);
end
