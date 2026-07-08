function [fig, fitStats] = plot_fz_vs_V(z_targets, Fz_over_m, SE_plot, ...
        d1, footArea, m_foot, rho_s, g_cm)
% PLOT_FZ_VS_V  Step 2 — F(z)/m vs V_hull(z)  (Kang modified Archimedes test).
%
%   Tests F(z)/m = slope * V_hull(z) forced through origin (stable depths only).
%   Theoretical Kang line removed — K_eff reported in annotation instead.
%   Vertical red dashed line marks stable shaft boundary in V space.
%
%   Symbol style matches plot_fz_vs_z:
%     Black edge, colored fill by regime (grey->mauve->crimson-pink)
%     SE bars at ALL depths (OLS-based SE passed in via SE_plot)
%
%   Inputs:
%       z_targets - [nZ x 1] depth targets [cm]
%       Fz_over_m - [nZ x 1] OLS intercepts (all depths)
%       SE_plot   - [nZ x 1] SE per depth (OLS-based, from plot_fz_vs_z fitStats)
%       d1        - shared d1 [cm] (annotation)
%       footArea  - struct from extract_foot_area_vs_depth
%       m_foot    - total dropped mass [g], default 60
%       rho_s     - bed BULK density [g/cm^3]; pass per-condition rho_bulk
%       g_cm      - gravity [cm/s^2], default 980

    if nargin < 6, m_foot = 60;   end
    if nargin < 7 || isempty(rho_s)
        warning('plot_fz_vs_V:noRho', ['rho_s (bed BULK density) not provided; falling ' ...
            'back to 1.52 g/cm^3. Pass per-condition rho_bulk from ' ...
            'get_substrate_properties(material,condition) for correct K_eff.']);
        rho_s = 1.52;
    end
    if nargin < 8, g_cm   = 980;  end

    z_targets = z_targets(:);
    Fz_over_m = Fz_over_m(:);
    SE_plot   = SE_plot(:);
    nZ        = numel(z_targets);
    z_min_fit = 1.30;
    z_detect  = 0.40;

    % Symbol colors — black edge, colored fill
    col_near   = [0.82 0.82 0.82];
    col_trans  = [0.55 0.35 0.45];
    col_stable = [0.55 0.05 0.25];
    ms = 9; lw_sym = 1.6;

    %% ── Compute V_hull(z) at all depths ──────────────────────────────────
    fa_z = footArea.depth_cm(:);
    fa_A = footArea.A_hull_sm(:);

    V_hull = nan(nZ, 1);
    for zi = 1:nZ
        z_q = z_targets(zi);
        if z_q > max(fa_z), continue; end
        z_fine     = linspace(0, z_q, 500);
        A_fine     = interp1(fa_z, fa_A, z_fine, 'pchip', 0);
        V_hull(zi) = trapz(z_fine, A_fine);
    end

    % V at stable boundary — for vertical transition line
    V_at_boundary = V_hull(find(z_targets >= z_min_fit, 1, 'first'));

    %% ── Fit: stable depths only ──────────────────────────────────────────
    is_stable = z_targets >= z_min_fit & isfinite(Fz_over_m) & ...
                Fz_over_m > 0 & isfinite(V_hull);
    V_fit = V_hull(is_stable);
    F_fit = Fz_over_m(is_stable);

    slope_fit = V_fit \ F_fit;
    pred_fit  = slope_fit .* V_fit;
    r2_fit    = 1 - sum((F_fit - pred_fit).^2) / ...
                    sum((F_fit - mean(F_fit)).^2);

    % K_eff and phi_eff
    K_phi_theory = 13.4;
    K_eff        = slope_fit * m_foot / (rho_s * g_cm);
    phi_test     = (20:0.1:50)' * pi/180;
    K_approx     = 2*(1+sin(phi_test))./(1-sin(phi_test)) .* exp(pi*tan(phi_test));
    [~, idx]     = min(abs(K_approx - K_eff));
    phi_eff_deg  = phi_test(idx) * 180/pi;

    fprintf('\n-- plot_fz_vs_V: Kang volume scaling ------------------------\n');
    fprintf('Fitted slope  = %.1f  R2=%.4f\n', slope_fit, r2_fit);
    fprintf('K_eff         = %.1f  (%.2fx theory K_phi=%.1f)\n', ...
        K_eff, K_eff/K_phi_theory, K_phi_theory);
    fprintf('phi_eff       = %.1f deg  (Delta=%.1f deg)\n', ...
        phi_eff_deg, phi_eff_deg-24);
    fprintf('-------------------------------------------------------------\n\n');

    %% ── Axis limits ──────────────────────────────────────────────────────
    V_max     = max(V_hull(isfinite(V_hull))) * 1.12;
    F_max     = max(F_fit + SE_plot(is_stable)) * 1.30;
    all_F     = Fz_over_m(isfinite(Fz_over_m) & isfinite(V_hull));
    y_lo      = min(0, min(all_F) * 1.15);

    %% ── Figure ───────────────────────────────────────────────────────────
    fig = figure('Name','Step 2: F(z)/m vs V_hull(z)', ...
                 'ToolBar','none','MenuBar','none');
    fig.Position = [160 100 760 600];
    ax = axes(fig, 'Position', [0.12 0.12 0.85 0.83]);
    hold(ax, 'on');

    % y=0 reference — black dashed
    yline(ax, 0, '--', 'Color',[0.10 0.10 0.10], 'LineWidth',1.8, ...
        'HandleVisibility','off');

    % Vertical transition line in V space — dark red dashed
    if isfinite(V_at_boundary)
        xline(ax, V_at_boundary, '--', 'Color',[0.65 0.05 0.05], ...
            'LineWidth', 2.5, 'HandleVisibility','off');
        text(ax, V_at_boundary - 0.03, F_max * 0.88, ...
            sprintf('geometric\ntransition'), ...
            'FontSize',18, 'Interpreter','latex', ...
            'Color',[0.65 0.05 0.05], ...
            'HorizontalAlignment','right', 'VerticalAlignment','top');
        text(ax, V_at_boundary + 0.03, F_max * 0.88, ...
            sprintf('stable shaft\n$V \\geq %.2f$ cm$^3$', V_at_boundary), ...
            'FontSize', 18, 'Interpreter', 'latex', ...
            'Color', [0.65 0.05 0.05], ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end

    % Fitted slope — dark blue dashed (Step 2 reference color)
    V_line = linspace(0, V_max, 300);
    plot(ax, V_line, slope_fit .* V_line, '--', ...
        'Color',[0.05 0.15 0.55], 'LineWidth', 2.5, ...
        'DisplayName', sprintf('F/m = %.0f V_hull,  R^2=%.2f', slope_fit, r2_fit));

    % Symbols — all depths, three regimes
    for zi = 1:nZ
        if ~isfinite(Fz_over_m(zi)) || ~isfinite(V_hull(zi)), continue; end
        z_q  = z_targets(zi);
        F_q  = Fz_over_m(zi);
        SE_q = SE_plot(zi);
        V_q  = V_hull(zi);

        if z_q < z_detect
            fc = col_near;
        elseif z_q < z_min_fit
            fc = col_trans;
        else
            fc = col_stable;
        end

        errorbar(ax, V_q, F_q, SE_q, SE_q, 0, 0, ...
            'o', 'Color','k', 'MarkerFaceColor', fc, ...
            'MarkerEdgeColor','k', 'MarkerSize', ms, ...
            'LineWidth', lw_sym, 'CapSize', 3, 'HandleVisibility','off');
    end

    set(ax, 'FontSize',13, 'Box','on', 'LineWidth',1.2, ...
        'XColor',[0 0 0], 'YColor',[0 0 0], ...
        'XMinorTick','on', 'YMinorTick','on', 'TickDir','in', ...
        'XLim',[0, V_max], 'YLim',[y_lo, F_max]);
    grid(ax,'off');
    xlabel(ax, '$V_{\rm hull}(z)$  (cm$^3$)', ...
        'FontSize',16, 'Interpreter','latex', 'Color',[0 0 0]);
    ylabel(ax, '$F(z)/m$  (cm s$^{-2}$)', ...
        'FontSize',16, 'Interpreter','latex', 'Color',[0 0 0]);



    %% ── Pack fitStats ────────────────────────────────────────────────────
    fitStats.slope_fit   = slope_fit;
    fitStats.r2_fit      = r2_fit;
    fitStats.K_eff       = K_eff;
    fitStats.K_phi_theory= K_phi_theory;
    fitStats.phi_eff_deg = phi_eff_deg;
    fitStats.V_hull      = V_hull;
    fitStats.z_min_fit   = z_min_fit;
    fitStats.V_at_boundary = V_at_boundary;
end