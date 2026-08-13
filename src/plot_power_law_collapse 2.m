function fig = plot_power_law_collapse(heights, cmap, d1, k_over_m)
% PLOT_POWER_LAW_COLLAPSE  Rigorous power-law proof for granular impact.
%
%   Main panel: v/v0 vs (1-z/z_max) on log-log axes.
%   Straight-line collapse proves v/v0 = (1-z/z_max)^alpha is a genuine
%   power law from continuous trajectory data.
%
%   Inset (bottom right): per-height alpha vs v0 on semi-log x-axis.
%   Bootstrap CIs, filled markers colored by height group (matching main
%   panel), log-linear trend alpha = a + b*ln(v0), four reference lines
%   matching main panel style exactly.
%
%   Inputs:
%       heights   - struct array from group_trials_by_height
%       cmap      - [nH x 3] colormap matching other figures
%       d1        - shared inertial length scale [cm]
%       k_over_m  - friction coefficient [s^-2]

    g_cm      = 980;
    dz        = 0.001;
    z_max_ode = 4.0;
    z_grid    = 0:dz:z_max_ode;
    nZ        = numel(z_grid);
    nH        = numel(heights);

    % Reference line colors — exact match to plot_normalized_collapse
    col_quad    = [0.05 0.45 0.45];     % teal dotted      — quadratic drag
    col_const   = [0.00 0.00 0.00];     % black dotted      — const drag
    col_lin     = [0.35 0.35 0.35];     % grey dash-dot     — linear resist
    col_fit     = [0.00 0.00 0.00];     % black dashed      — global fit
    col_inertial= [0.05 0.15 0.55];     % dark blue dashed  — pure inertial (main only)
    col_kats    = [0.65 0.05 0.05];     % dark red dashed   — Katsuragi sphere

    % Separate darkest blues in cmap (match plot_normalized_collapse)
    cmap(1,:) = [0.10 0.10 0.55];
    cmap(2,:) = [0.35 0.55 0.95];

    %% ── Collect trajectories + per-height alpha ──────────────────────────
    all_xi   = [];
    all_vi   = [];
    trial_xi = {};
    trial_vi = {};
    trial_col= [];
    ti       = 0;

    % Per-height bootstrap storage
    alpha_j   = nan(nH, 1);
    alpha_jlo = nan(nH, 1);
    alpha_jhi = nan(nH, 1);
    v0vals    = nan(nH, 1);

    for j = 1:nH
        v0vals(j) = heights(j).v0_mean;
        zj = []; vj = []; gj_z = {}; gj_v = {}; gc = 0;

        for i = 1:heights(j).nTrials
            k_t    = heights(j).trials(i).kinematics;
            s_t    = heights(j).trials(i).scalars;
            idx    = k_t.impact_index:k_t.stopFrame;
            z_meas = k_t.z_smooth(idx);
            v_meas = k_t.v_smooth(idx);
            v0_t   = s_t.v0_cm_s;
            zmax_t = s_t.d_final_cm;

            if zmax_t <= 0 || ~isfinite(zmax_t) || v0_t <= 0, continue; end

            % Main panel: 1-z/zmax vs v/v0
            xi = 1 - z_meas ./ zmax_t;
            vi = v_meas ./ v0_t;
            ok = xi > 0.01 & xi <= 1.0 & vi > 0.005 & vi <= 1.0 & ...
                 isfinite(xi) & isfinite(vi);
            if sum(ok) < 5, continue; end

            all_xi = [all_xi; xi(ok)];
            all_vi = [all_vi; vi(ok)];
            ti = ti + 1;
            trial_xi{ti}  = xi(ok);
            trial_vi{ti}  = vi(ok);
            trial_col(ti) = j;

            % Inset: z/zmax basis for per-height alpha
            z_n   = z_meas ./ zmax_t;
            v_n   = v_meas ./ v0_t;
            valid = isfinite(z_n) & isfinite(v_n) & z_n >= 0 & z_n < 0.99 & v_n > 0;
            if ~any(valid), continue; end
            zj = [zj; z_n(valid)];
            vj = [vj; v_n(valid)];
            gc = gc + 1;
            gj_z{gc} = z_n(valid);
            gj_v{gc} = v_n(valid);
        end

        if numel(zj) < 10, continue; end
        alpha_j(j) = log(1 - zj) \ log(vj);

        % Bootstrap CI per height group
        if gc >= 2
            ab = nan(300,1);
            for kb = 1:300
                ib  = randi(gc, gc, 1);
                zb  = vertcat(gj_z{ib});
                vb  = vertcat(gj_v{ib});
                mb  = zb < 0.99 & vb > 0;
                if sum(mb) < 5, continue; end
                ab(kb) = log(1 - zb(mb)) \ log(vb(mb));
            end
            ab = ab(isfinite(ab));
            if numel(ab) > 10
                alpha_jlo(j) = prctile(ab, 5);
                alpha_jhi(j) = prctile(ab, 95);
            end
        end
    end

    %% ── Global power law fit (log-log space, forced through origin) ──────
    log_xi    = log(all_xi);
    log_vi    = log(all_vi);
    alpha_fit = log_xi \ log_vi;

    xi_line = logspace(log10(0.008), log10(1.0), 400);
    r2_pow  = 1 - sum((log_vi - alpha_fit.*log_xi).^2) / ...
                  sum((log_vi - mean(log_vi)).^2);

    %% ── Log-linear inset fit: alpha_j = a + b*ln(v0) ────────────────────
    ok        = isfinite(alpha_j) & isfinite(v0vals);
    X_ln      = [ones(sum(ok),1), log(v0vals(ok))];
    p_ln      = X_ln \ alpha_j(ok);
    v0ref     = linspace(min(v0vals(ok))*0.92, max(v0vals(ok))*1.08, 200);
    alpha_trend = p_ln(1) + p_ln(2) .* log(v0ref);

    alpha_pred = p_ln(1) + p_ln(2) .* log(v0vals(ok));
    ss_res     = sum((alpha_j(ok) - alpha_pred).^2);
    ss_tot     = sum((alpha_j(ok) - mean(alpha_j(ok))).^2);
    r2_loglin  = 1 - ss_res / ss_tot;

    fprintf('\n-- plot_power_law_collapse ----------------------------------\n');
    fprintf('Global alpha = %.4f  (R2 = %.4f)\n', alpha_fit, r2_pow);
    fprintf('Log-linear fit: alpha = %.3f + %.3f*ln(v0)  R2 = %.3f\n', ...
        p_ln(1), p_ln(2), r2_loglin);
    fprintf('Foot sits %.1f%% between pure-inertial and Katsuragi limits\n', ...
        100*(alpha_fit-0.5)/(2/3-0.5));
    fprintf('------------------------------------------------------------\n\n');

    %% ── ODE predictions per height group ─────────────────────────────────
    ode_xi = cell(nH,1);
    ode_vi = cell(nH,1);
    for j = 1:nH
        v0_j      = heights(j).v0_mean;
        v2_ode    = zeros(1,nZ);
        v2_ode(1) = v0_j^2;
        for zi = 1:nZ-1
            z_i  = z_grid(zi);
            v2_i = v2_ode(zi);
            dv2  = 2*(g_cm - k_over_m*z_i - v2_i/d1);
            v2_ode(zi+1) = v2_i + dv2*dz;
            if v2_ode(zi+1) <= 0
                v2_ode(zi+1:end) = 0; break;
            end
        end
        v_ode    = sqrt(max(v2_ode,0));
        idx_stop = find(v2_ode <= 0, 1, 'first');
        d_stop   = heights(j).d_mean;
        if ~isempty(idx_stop), d_stop = z_grid(idx_stop); end
        if d_stop <= 0, continue; end
        xi_o = 1 - z_grid./d_stop;
        vi_o = v_ode./v0_j;
        ok_o = xi_o > 0.008 & xi_o <= 1.0 & vi_o > 0.005;
        ode_xi{j} = xi_o(ok_o);
        ode_vi{j} = vi_o(ok_o);
    end

    %% ── Figure ───────────────────────────────────────────────────────────
    fig = figure('Name','Power-law collapse', ...
                 'ToolBar','none','MenuBar','none');
    fig.Position = [120 80 920 640];

    ax = axes(fig, 'Position',[0.10 0.12 0.86 0.83]);
    hold(ax,'on');
    set(ax, 'XScale','log','YScale','log', ...
        'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in', ...
        'XLim',[0.008 1.2],'YLim',[0.008 1.2]);

    % Measured trajectories — colored by height group
    for t = 1:numel(trial_xi)
        plot(ax, trial_xi{t}, trial_vi{t}, '-', ...
            'Color',[cmap(trial_col(t),:) 0.18], 'LineWidth',0.9, ...
            'HandleVisibility','off');
    end

    % ODE predictions — dotted, same color scheme
    for j = 1:nH
        if isempty(ode_xi{j}), continue; end
        plot(ax, ode_xi{j}, ode_vi{j}, ':', ...
            'Color',[cmap(j,:) 0.65], 'LineWidth',1.5, ...
            'HandleVisibility','off');
    end

    % Reference lines — same style as plot_normalized_collapse
    plot(ax, xi_line, xi_line.^1.0, ':', ...
        'Color',col_quad,'LineWidth',2.4, ...
        'DisplayName','$\alpha = 1.00$ (quadratic drag)');
    plot(ax, xi_line, xi_line.^0.5, ':', ...
        'Color',col_const,'LineWidth',2.8, ...
        'DisplayName','$\alpha = 0.50$ (constant drag)');
    plot(ax, xi_line, xi_line.^(2/3), '-.', ...
        'Color',col_lin,'LineWidth',2.8, ...
        'DisplayName','$\alpha = 0.67$ (Katsuragi sphere)');
    plot(ax, xi_line, xi_line.^alpha_fit, '--', ...
        'Color',col_fit,'LineWidth',3.4, ...
        'DisplayName',sprintf('$\\alpha = %.3f$ (jerboa foot fit)', alpha_fit));

    xlabel(ax, '$1 - z/z_{\max}$','FontSize',22, ...
        'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax, '$v/v_0$','FontSize',22, ...
        'Interpreter','latex','Color',[0 0 0]);
    legend(ax,'show','FontSize',10,'Interpreter','latex','Box','on', ...
        'EdgeColor',[0.25 0.25 0.25],'Location','southeast');


    %% ── Inset: alpha_j vs v0, semi-log x, bottom right ──────────────────
    ax_in = axes(fig, 'Position',[0.62 0.18 0.28 0.32]);
    hold(ax_in,'on');

    % Reference horizontals — exact match to main panel styles
    yline(ax_in, 1.00, ':', 'Color',col_quad,  'LineWidth',2.8, ...
        'HandleVisibility','off');
    yline(ax_in, 0.50, ':', 'Color',col_const, 'LineWidth',3.2, ...
        'HandleVisibility','off');
    yline(ax_in, 2/3,  '-.','Color',col_lin,   'LineWidth',3.2, ...
        'HandleVisibility','off');
    yline(ax_in, alpha_fit,'--','Color',col_fit,'LineWidth',3.4, ...
        'HandleVisibility','off');

    % Log-linear trend line
    plot(ax_in, v0ref, alpha_trend, 'k-','LineWidth',2.5, ...
        'HandleVisibility','off');

    % Per-height markers with bootstrap CI — colored fill, black edge
    ok_idx = find(ok);
    for jj = 1:numel(ok_idx)
        j = ok_idx(jj);
        [~, marker] = get_height_style(heights(j).h_cm);
        lo = alpha_j(j) - alpha_jlo(j);
        hi = alpha_jhi(j) - alpha_j(j);
        if ~isfinite(lo), lo = 0; end
        if ~isfinite(hi), hi = 0; end
        errorbar(ax_in, v0vals(j), alpha_j(j), lo, hi, ...
            marker, ...
            'Color',           'k', ...
            'MarkerFaceColor', cmap(j,:), ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize',      10, ...
            'LineWidth',       1.2, ...
            'CapSize',         3, ...
            'HandleVisibility','off');
    end

    % R² annotation
    text(ax_in, 0.97, 0.08, sprintf('$R^2 = %.3f$', r2_loglin), ...
        'Units','normalized','Interpreter','latex','FontSize',9, ...
        'HorizontalAlignment','right','Color',[0.2 0.2 0.2]);

    alpha_jhi_max = max(alpha_jhi(ok & isfinite(alpha_jhi)));
    if isempty(alpha_jhi_max) || ~isfinite(alpha_jhi_max)
        alpha_jhi_max = 1.0;
    end
    set(ax_in, 'XScale','log', ...
        'FontSize',10,'Box','on','LineWidth',1.0, ...
        'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','in', ...
        'XMinorTick','on','YMinorTick','on', ...
        'XLim',[min(v0vals(ok))*0.92, max(v0vals(ok))*1.08], ...
        'YLim',[0, max(1.05, alpha_jhi_max*1.05)]);
    grid(ax_in,'off');
    xlabel(ax_in,'$v_0$  (cm s$^{-1}$)','FontSize',18, ...
        'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax_in,'Fitted $\alpha$','FontSize',18, ...
        'Interpreter','latex','Color',[0 0 0]);
end