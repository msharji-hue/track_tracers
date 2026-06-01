function [fig, d0] = plot_depth_vs_v0(heights, cmap, footArea, d1, k_over_m)
% PLOT_DEPTH_VS_V0  Final penetration depth vs impact speed.
%   Replicates Katsuragi & Durian 2007 Fig. 2b.
%   Always fits and returns d0 from the (d0²H)^(1/3) depth-scaling law.
%   Forward model overlay shown only when footArea, d1, k_over_m supplied.
%
%   Outputs:
%     fig - figure handle
%     d0  - depth scale from (d0²H)^(1/3) fit [cm] — feeds plot_fz_vs_z

    doForward = nargin >= 3 && ~isempty(footArea);
    g_cm      = 980;

    % ── Collect per-trial scalars ─────────────────────────────────────────
    v0_all = [];
    d_all  = [];
    H_all  = [];

    for j = 1:numel(heights)
        for i = 1:heights(j).nTrials
            s  = heights(j).trials(i).scalars;
            h  = heights(j).h_cm;
            d  = s.d_final_cm;
            v0_all(end+1) = s.v0_cm_s;
            d_all(end+1)  = d;
            H_all(end+1)  = h + d;   % total drop distance H = h + d
        end
    end

    ss_tot = sum((d_all - mean(d_all)).^2);

    % ── Linear fit: d = d0_lin + alpha*|v0| ──────────────────────────────
    X         = [ones(numel(v0_all),1), v0_all(:)];
    coeffs    = X \ d_all(:);
    d0_lin    = coeffs(1);
    alpha_fit = coeffs(2);
    r2_lin    = 1 - sum((d_all - (d0_lin + alpha_fit.*v0_all)).^2) / ss_tot;

    % ── Depth-scaling fit: d = (d0² * H)^(1/3) ───────────────────────────
    % Rearranged: d0² = d³ / H  →  d0 = sqrt(mean(d³ ./ H))
    d0_sq = mean(d_all.^3 ./ H_all);
    d0    = sqrt(d0_sq);

    d_scaling_pred = (d0^2 .* H_all).^(1/3);
    r2_scaling     = 1 - sum((d_all - d_scaling_pred).^2) / ss_tot;

    % ── Katsuragi velocity-based model: d = (d1/2)*ln(1 + v0²/(g*d1)) ───
    kats_model = @(d1k, v0) (d1k/2) .* log(1 + v0.^2 ./ (g_cm * d1k));
    d1_kats    = fminsearch(@(d1k) sum((d_all - kats_model(d1k, v0_all)).^2), 5);
    r2_kats    = 1 - sum((d_all - kats_model(d1_kats, v0_all)).^2) / ss_tot;

    fprintf('\n-- d vs v0 --------------------------------------------------\n');
    fprintf('Linear:         d0=%.3f cm  alpha=%.5f s/cm  R2=%.4f\n', ...
        d0_lin, alpha_fit, r2_lin);
    fprintf('Depth-scaling:  d0=%.3f cm  (d0²H)^1/3       R2=%.4f\n', ...
        d0, r2_scaling);
    fprintf('Katsuragi vel:  d1=%.3f cm                    R2=%.4f\n', ...
        d1_kats, r2_kats);
    fprintf('-------------------------------------------------------------\n\n');

    % ── Forward model (optional) ──────────────────────────────────────────
    gammas      = [0, 0.3, 0.5, 1.0];
    gamma_cols  = [0.20 0.55 0.85;
                   0.15 0.70 0.35;
                   0.85 0.50 0.05;
                   0.75 0.10 0.10];
    gamma_lines = {'-.', '-', '--', ':'};

    v0_line = linspace(min(v0_all)*0.90, max(v0_all)*1.05, 60);
    d_fwd   = nan(numel(gammas), numel(v0_line));
    r2_fwd  = nan(numel(gammas), 1);

    if doForward
        max_d_exp = max(d_all) * 1.1;
        in_range  = footArea.depth_cm <= max_d_exp;
        A_ref     = mean(footArea.A_bare_sm(in_range), 'omitnan');

        fprintf('-- Forward model --------------------------------------------\n');
        fprintf('A_ref = %.4f cm²\nd1 = %.3f cm\nk/m = %.1f s^-2\n\n', ...
            A_ref, d1, k_over_m);

        for q = 1:numel(gammas)
            A_eff_q = footArea.A_bare_sm + gammas(q) .* footArea.A_gap_sm;
            A_eff_q = max(A_eff_q, 0);
            for vi = 1:numel(v0_line)
                d_fwd(q,vi) = run_forward_model(v0_line(vi), d1, k_over_m, ...
                                                A_eff_q, footArea.depth_cm, A_ref);
            end
            d_pred_trials = arrayfun(@(v) run_forward_model(v, d1, k_over_m, ...
                                A_eff_q, footArea.depth_cm, A_ref), v0_all);
            r2_fwd(q) = 1 - sum((d_all - d_pred_trials).^2) / ss_tot;
            fprintf('gamma=%.1f:  R2=%.4f\n', gammas(q), r2_fwd(q));
        end
        fprintf('-------------------------------------------------------------\n\n');
    end

    % ── Reference curves ──────────────────────────────────────────────────
    H_line    = v0_line.^2 ./ (2*g_cm) + (d0_lin + alpha_fit.*v0_line);
    d_scale   = (d0^2 .* H_line).^(1/3);
    d_lin_ref = d0_lin + alpha_fit .* v0_line;
    d_kats    = kats_model(d1_kats, v0_line);

    % ── Figure ────────────────────────────────────────────────────────────
    fig = figure('Name','d vs v0','ToolBar','none','MenuBar','none');
    fig.Position = [100 100 660 520];
    ax = axes(fig, 'Position', [0.12 0.13 0.60 0.82]);
    hold(ax,'on');

    % Reference curves
plot(ax, v0_line, d_lin_ref, '--', 'Color', [0.25 0.25 0.25], ...
    'LineWidth', 2.0, 'HandleVisibility', 'off');

% Depth-scaling curve intentionally not displayed
% d0 is still computed and returned for downstream use

plot(ax, v0_line, d_kats, ':', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.6, 'HandleVisibility', 'off');

    % Forward model curves
    if doForward
        for q = 1:numel(gammas)
            plot(ax, v0_line, d_fwd(q,:), gamma_lines{q}, ...
                'Color', gamma_cols(q,:), 'LineWidth', 2.0, ...
                'HandleVisibility', 'off');
        end
    end

    % Data points per height group
    for j = 1:numel(heights)
        hg          = heights(j);
        [~, marker] = get_height_style(hg.h_cm);
        errorbar(ax, hg.v0_mean, hg.d_mean, hg.d_std, hg.d_std, ...
            hg.v0_std, hg.v0_std, marker, ...
            'Color', cmap(j,:), 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', cmap(j,:), 'MarkerSize', 9, 'LineWidth', 1.8, ...
            'HandleVisibility', 'off');
    end

    set(ax,'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in', ...
        'XLim',[70, max(v0_all)*1.05], 'YLim',[0, max(d_all)*1.15]);
    grid(ax,'off');
    xlabel(ax,'$v_0$  (cm s$^{-1}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax,'$d$  (cm)',            'FontSize',16,'Interpreter','latex','Color',[0 0 0]);

    % ── Legend (dummy axes) ───────────────────────────────────────────────
    ax2 = axes(fig, 'Position', ax.Position, 'Visible','off');
    hold(ax2,'on');
    
    plot(ax2, nan, nan, '--', 'Color',[0.25 0.25 0.25], 'LineWidth',2.0, ...
        'DisplayName', sprintf('$d = %.2f + %.4f v_0$', d0_lin, alpha_fit));
    
    plot(ax2, nan, nan, ':', 'Color',[0.35 0.35 0.35], 'LineWidth',1.6, ...
        'DisplayName', sprintf('Katsuragi velocity fit, $R^2=%.3f$', r2_kats));
    
    if doForward
        for q = 1:numel(gammas)
            plot(ax2, nan, nan, gamma_lines{q}, ...
                'Color', gamma_cols(q,:), 'LineWidth',1.8, ...
                'DisplayName', sprintf('$\\gamma=%.1f$, $R^2=%.3f$', ...
                    gammas(q), r2_fwd(q)));
        end
    end
    
    lgd2 = legend(ax2,'show', ...
        'FontSize',8.5, ...
        'Box','on', ...
        'Interpreter','latex', ...
        'EdgeColor',[0.15 0.15 0.15], ...
        'Color',[1 1 1]);
    
    lgd2.Location = 'none';
    
    if doForward
        lgd2.Position = [0.735 0.58 0.245 0.25];
    else
        lgd2.Position = [0.57 0.16 0.34 0.10];
    end
    % ── Annotation ────────────────────────────────────────────────────────
    text(ax, 0.03, 0.97, ...
        sprintf(['$R^2_{\\rm lin} = %.3f$\n' ...
                 '$R^2_{\\rm K.vel} = %.3f$\n' ...
                 '$d_0^{\\rm scale} = %.2f$ cm'], ...
            r2_lin, r2_kats, d0), ...
        'FontSize',9,'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'Interpreter','latex','FontAngle','italic', ...
        'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25], ...
        'LineWidth',1.0,'Margin',4);
    end