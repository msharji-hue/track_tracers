function fig = plot_depth_vs_v0(heights, cmap)
% PLOT_DEPTH_VS_V0  Replicates Katsuragi & Durian 2007 Fig. 2b.
%   Two curves only — matching paper exactly:
%     grey dashed : d = d0 + alpha*|v0|           [de Bruyn & Walsh 2004]
%     black dotted: d = (d1/2)*ln(1+v0^2/(g*d1)) [Katsuragi & Durian 2007]

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
    kats_model = @(d1, v0) (d1/2) .* log(1 + v0.^2 ./ (g_cm * d1));
    d1_fit     = fminsearch(@(d1) sum((d_all - kats_model(d1, v0_all)).^2), 5);
    r2_kats    = 1 - sum((d_all - kats_model(d1_fit, v0_all)).^2) / ss_tot;

    fprintf('\n── d vs v0 fit ──────────────────────────────────\n');
    fprintf('Linear:    d0=%.3f cm  alpha=%.5f s/cm  R2=%.4f\n', d0_fit, alpha_fit, r2_lin);
    fprintf('Katsuragi: d1=%.3f cm                    R2=%.4f\n', d1_fit, r2_kats);
    fprintf('─────────────────────────────────────────────────\n\n');

    % ── Reference curves ──────────────────────────────────────────────────
    v0_min  = max(0, min(v0_all)*0.85);
    v0_line = linspace(v0_min, max(v0_all)*1.05, 200);
    d_lin   = d0_fit + alpha_fit .* v0_line;
    d_kats  = kats_model(d1_fit, v0_line);

    % ── Figure ────────────────────────────────────────────────────────────
    fig = figure('Name','d vs v0','ToolBar','none','MenuBar','none');
    fig.Position = [100 100 620 520];
    ax = axes(fig, 'Position', [0.12 0.13 0.62 0.82]);
    hold(ax, 'on');

    plot(ax, v0_line, d_lin, '--', 'Color', [0.55 0.55 0.55], ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');
    plot(ax, v0_line, d_kats, ':', 'Color', [0 0 0], ...
        'LineWidth', 2.2, 'HandleVisibility', 'off');

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
        'XLim',[v0_min, max(v0_all)*1.05]);
    grid(ax,'off');
    xlabel(ax,'$v_0$  (cm s$^{-1}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax,'$d$  (cm)',            'FontSize',16,'Interpreter','latex','Color',[0 0 0]);

    % Model legend
    ax2 = axes(fig, 'Position', ax.Position, 'Visible', 'off');
    hold(ax2, 'on');
    plot(ax2, nan, nan, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.8, ...
        'DisplayName', 'd_0 + \alpha|v_0|');
    plot(ax2, nan, nan, ':',  'Color', [0 0 0],           'LineWidth', 2.2, ...
        'DisplayName', 'Model  [Katsuragi and Durian, 2007]');
    lgd2 = legend(ax2,'show','FontSize',11,'Box','on','Interpreter','tex', ...
        'EdgeColor',[0.15 0.15 0.15],'LineWidth',1.2,'Color',[1 1 1]);
    lgd2.Location = 'none';
    lgd2.Position = [0.14 0.76 0.50 0.13];

    text(ax, 0.03, 0.97, ...
        sprintf(['d_0 = %.2f cm,  \\alpha = %.4f s cm^{-1}\n' ...
                 'd_1 = %.2f cm  (depth-scaling fit)\n'       ...
                 'R^2_{lin}  = %.4f\n'                         ...
                 'R^2_{kats} = %.4f'], ...
                d0_fit, alpha_fit, d1_fit, r2_lin, r2_kats), ...
        'FontSize',10,'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'Color',[0.10 0.10 0.10],'FontAngle','italic','Interpreter','tex', ...
        'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25], ...
        'LineWidth',1.0,'Margin',4);
end
