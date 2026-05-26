function fig = plot_depth_vs_v0(heights, cmap, footArea, d1, k_over_m)
% PLOT_DEPTH_VS_V0  Replicates Katsuragi & Durian 2007 Fig. 2b.
%   Two reference curves (linear + Katsuragi) plus optional forward model
%   overlay for multiple gamma values if footArea struct is provided.
%
%   Inputs:
%     heights    - struct array from group_trials_by_height
%     cmap       - colormap
%     footArea   - (optional) struct from extract_foot_area_vs_depth
%     d1         - (optional) kinematic d1 (cm), required if footArea given
%     k_over_m   - (optional) k/m (s^-2) from plot_fz_vs_z, required if footArea given

    doForward = nargin >= 3 && ~isempty(footArea);

    g_cm = 980;

    % ── Collect all individual trial scalars ──────────────────────────────
    v0_all = [];
    d_all  = [];
    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials
            v0_all(end+1) = heights(j).trials(i).scalars.v0_cm_s;
            d_all(end+1)  = heights(j).trials(i).scalars.d_final_cm;
        end
    end
    ss_tot = sum((d_all - mean(d_all)).^2);

    % ── Linear fit: d = d0 + alpha*|v0| ──────────────────────────────────
    X         = [ones(numel(v0_all),1), v0_all(:)];
    coeffs    = X \ d_all(:);
    d0_fit    = coeffs(1);
    alpha_fit = coeffs(2);
    r2_lin    = 1 - sum((d_all-(d0_fit+alpha_fit.*v0_all)).^2) / ss_tot;

    % ── Katsuragi model: d = (d1/2)*ln(1 + v0^2/(g*d1)) ─────────────────
    kats_model = @(d1k, v0) (d1k/2) .* log(1 + v0.^2 ./ (g_cm * d1k));
    d1_kats    = fminsearch(@(d1k) sum((d_all - kats_model(d1k, v0_all)).^2), 5);
    r2_kats    = 1 - sum((d_all - kats_model(d1_kats, v0_all)).^2) / ss_tot;

    fprintf('\n── d vs v0 fit ──────────────────────────────────\n');
    fprintf('Linear:    d0=%.3f cm  alpha=%.5f s/cm  R2=%.4f\n', d0_fit, alpha_fit, r2_lin);
    fprintf('Katsuragi: d1=%.3f cm                    R2=%.4f\n', d1_kats, r2_kats);

    % ── Forward model predictions for multiple gammas ─────────────────────
    gammas      = [0, 0.3, 0.5, 1.0];
    gamma_cols  = [0.20 0.55 0.85;   % blue   — γ=0
                   0.15 0.70 0.35;   % green  — γ=0.3
                   0.85 0.50 0.05;   % orange — γ=0.5
                   0.75 0.10 0.10];  % red    — γ=1.0
    gamma_lines = {'-.',  '-',  '--',  ':'};

    v0_line = linspace(min(v0_all)*0.90, max(v0_all)*1.05, 60);
    d_fwd   = nan(numel(gammas), numel(v0_line));
    r2_fwd  = nan(numel(gammas), 1);

    if doForward
        % Reference area: mean A_eff over the experimental depth range
        max_d_exp = max(d_all) * 1.1;
        in_range  = footArea.depth_cm <= max_d_exp;
        A_ref     = mean(footArea.A_bare_sm(in_range), 'omitnan');

        fprintf('\n── Forward model ────────────────────────────────\n');
        fprintf('A_ref = %.4f cm²  (mean bare area, z < %.2f cm)\n', A_ref, max_d_exp);
        fprintf('d1    = %.3f cm\n', d1);
        fprintf('k/m   = %.1f s^-2\n\n', k_over_m);

        for q = 1:numel(gammas)
            g_val    = gammas(q);
            A_eff_q  = footArea.A_bare_sm + g_val .* footArea.A_gap_sm;
            A_eff_q  = max(A_eff_q, 0);

            for vi = 1:numel(v0_line)
                d_fwd(q,vi) = run_forward_model(v0_line(vi), d1, k_over_m, ...
                                                A_eff_q, footArea.depth_cm, A_ref);
            end

            % R² against per-trial data
            d_pred_trials = arrayfun(@(v) run_forward_model(v, d1, k_over_m, ...
                                A_eff_q, footArea.depth_cm, A_ref), v0_all);
            r2_fwd(q) = 1 - sum((d_all - d_pred_trials).^2) / ss_tot;

            fprintf('gamma=%.1f:  R2=%.4f\n', g_val, r2_fwd(q));
        end
        fprintf('─────────────────────────────────────────────────\n\n');
    end

    % ── Reference curves ──────────────────────────────────────────────────
    d_lin  = d0_fit + alpha_fit .* v0_line;
    d_kats = kats_model(d1_kats, v0_line);

    % ── Figure ────────────────────────────────────────────────────────────
    fig = figure('Name','d vs v0','ToolBar','none','MenuBar','none');
    fig.Position = [100 100 660 520];
    ax = axes(fig, 'Position', [0.12 0.13 0.60 0.82]);
    hold(ax, 'on');

    % Reference curves
    plot(ax, v0_line, d_lin, '--', 'Color', [0.55 0.55 0.55], ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');
    plot(ax, v0_line, d_kats, ':', 'Color', [0 0 0], ...
        'LineWidth', 2.2, 'HandleVisibility', 'off');

    % Forward model curves
    if doForward
        for q = 1:numel(gammas)
            plot(ax, v0_line, d_fwd(q,:), gamma_lines{q}, ...
                'Color', gamma_cols(q,:), 'LineWidth', 2.0, ...
                'HandleVisibility', 'off');
        end
    end

    % Data points
    for j = 1:numel(heights)
        hg          = heights(j);
        col         = cmap(j,:);
        [~, marker] = get_height_style(hg.h_cm);
        errorbar(ax, hg.v0_mean, hg.d_mean, hg.d_std, hg.d_std, ...
            hg.v0_std, hg.v0_std, marker, ...
            'Color', col, 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', col, 'MarkerSize', 9, 'LineWidth', 1.8, ...
            'HandleVisibility', 'off');
    end

    set(ax, 'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in', ...
        'XLim',[0, max(v0_all)*1.05], 'YLim',[0, max(d_all)*1.15]);
    grid(ax,'off');
    xlabel(ax,'$v_0$  (cm s$^{-1}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax,'$d$  (cm)',            'FontSize',16,'Interpreter','latex','Color',[0 0 0]);

    % ── Legend (dummy axes, right panel) ──────────────────────────────────
    ax2 = axes(fig, 'Position', ax.Position, 'Visible', 'off');
    hold(ax2, 'on');

    % Reference lines
    plot(ax2, nan, nan, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.8, ...
        'DisplayName', '$d_0 + \alpha|v_0|$');
    plot(ax2, nan, nan, ':',  'Color', [0 0 0], 'LineWidth', 2.2, ...
        'DisplayName', 'Katsuragi \& Durian (2007)');

    % Forward model lines
    if doForward
        for q = 1:numel(gammas)
            plot(ax2, nan, nan, gamma_lines{q}, ...
                'Color', gamma_cols(q,:), 'LineWidth', 2.0, ...
                'DisplayName', sprintf('Forward model,  $\\gamma = %.1f$  ($R^2=%.3f$)', ...
                    gammas(q), r2_fwd(q)));
        end
    end

    lgd2 = legend(ax2,'show','FontSize',10,'Box','on','Interpreter','latex', ...
        'EdgeColor',[0.15 0.15 0.15],'LineWidth',1.0,'Color',[1 1 1]);
    lgd2.Location = 'none';
    if doForward
        lgd2.Position = [0.74 0.55 0.24 0.35];
    else
        lgd2.Position = [0.14 0.76 0.50 0.13];
    end

    % ── Stats annotation ──────────────────────────────────────────────────
    txt = sprintf(['$d_0 = %.2f$ cm,  $\\alpha = %.4f$ s cm$^{-1}$\n' ...
                   '$d_1 = %.2f$ cm  (depth-scaling fit)\n'            ...
                   '$R^2_{\\rm lin}  = %.4f$\n'                         ...
                   '$R^2_{\\rm kats} = %.4f$'], ...
                  d0_fit, alpha_fit, d1_kats, r2_lin, r2_kats);
    text(ax, 0.03, 0.97, txt, ...
        'FontSize',10,'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'Color',[0.10 0.10 0.10],'FontAngle','italic','Interpreter','latex', ...
        'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25], ...
        'LineWidth',1.0,'Margin',4);
end