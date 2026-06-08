function [fig, fitStats] = plot_fz_vs_Az(z_targets, Fz_over_m, Fz_pooled_se, ...
        d1, k_over_m, footArea, m_foot, rho_s, g_cm)
% PLOT_FZ_VS_AZ  Step 3 — F(z)/m vs A_hull(z)*z (LFFM-inspired geometry).
%
%   Main panel: stable-shaft OLS intercepts vs A_hull(z)*z.
%   Two reference lines for your foot (dark green solid/dashed).
%   Two reference lines for comparison (black dotted = Katsuragi dynamic,
%   grey dotted = Kang quasistatic theory).
%
%   Inset (bottom right): k/m vs K_eff for each known geometry.
%   Shows how the jerboa foot occupies a different point in (k/m, K_eff)
%   space compared to the Katsuragi sphere and Kang theory. Future foot
%   geometries (different toe spacings) will add new points here.
%
%   Reference line colors (consistent across Step figures):
%     Dark green solid  [0.05 0.40 0.15] — F ~ C*Az       (your foot, linear)
%     Dark green dashed [0.05 0.40 0.15] — F ~ C*[Az]^p   (your foot, power law)
%     Black dotted      [0.05 0.05 0.05] — Katsuragi sphere (Approach B, dynamic)
%     Grey dotted       [0.55 0.55 0.55] — Kang theory (Approach A, quasistatic)
%
%   Inputs:
%       z_targets    - [nZ x 1] depth targets [cm]
%       Fz_over_m    - [nZ x 1] OLS intercepts (all depths)
%       Fz_pooled_se - [nZ x 1] SE per depth (from fitStats7.SE_plot)
%       d1           - shared d1 [cm]
%       k_over_m     - friction coefficient [s^-2] from Step 1
%       footArea     - struct with depth_cm, A_hull_sm
%       m_foot       - total dropped mass [g], default 60
%       rho_s        - packing density [g/cm^3], default 1.52
%       g_cm         - gravity [cm/s^2], default 980

    if nargin < 7, m_foot = 60;   end
    if nargin < 8, rho_s  = 1.52; end
    if nargin < 9, g_cm   = 980;  end

    z_targets    = z_targets(:);
    Fz_over_m    = Fz_over_m(:);
    Fz_pooled_se = Fz_pooled_se(:);
    nZ           = numel(z_targets);
    z_min_fit    = 1.30;

    % ── Katsuragi & Durian (2007) reference parameters ────────────────────
    m_sphere    = 28.4;                         % g
    D_sphere    = 2.54;                         % cm
    A_sphere    = pi*(D_sphere/2)^2;            % 5.067 cm²
    km_sphere   = 1200;                         % s^-2 (from their Fig 3b)
    d1_sphere   = 8.7;                          % cm
    slope_B     = km_sphere / A_sphere;         % 236.8 — Approach B (dynamic)
    K_eff_sphere= slope_B * m_sphere / (rho_s * g_cm);  % 4.52

    % ── Kang theory reference (Approach A) ───────────────────────────────
    K_phi_theory = 13.4;
    slope_A      = K_phi_theory * rho_s * g_cm / m_foot;  % 332.7

    % Symbol colors
    col_stable = [0.55 0.05 0.25];   % dark crimson-pink
    col_green  = [0.05 0.40 0.15];   % dark green — your foot fits
    col_kats   = [0.05 0.05 0.05];   % black dotted — Katsuragi sphere
    col_kang   = [0.55 0.55 0.55];   % grey dotted — Kang theory
    ms = 9; lw_sym = 1.6;

    %% ── Compute A_hull and Az at stable depths ───────────────────────────
    fa_z = footArea.depth_cm(:);
    fa_A = footArea.A_hull_sm(:);

    A_hull = nan(nZ, 1);
    Az     = nan(nZ, 1);
    for zi = 1:nZ
        z_q = z_targets(zi);
        if z_q > max(fa_z), continue; end
        A_hull(zi) = interp1(fa_z, fa_A, z_q, 'pchip', 0);
        Az(zi)     = A_hull(zi) * z_q;
    end

    % Stable depth arrays
    is_stable = z_targets >= z_min_fit & ...
                isfinite(Fz_over_m) & Fz_over_m > 0 & ...
                isfinite(Az) & Az > 0;
    Az_fit = Az(is_stable);
    F_fit  = Fz_over_m(is_stable);
    SE_fit = Fz_pooled_se(is_stable);

    if numel(Az_fit) < 2
        error('plot_fz_vs_Az: fewer than 2 stable-depth points available.');
    end

    fprintf('\nStable depths in fit: %d points\n', numel(Az_fit));
    fprintf('  Az range: %.4f - %.4f cm^3\n', min(Az_fit), max(Az_fit));
    fprintf('  F  range: %.1f - %.1f cm/s^2\n', min(F_fit), max(F_fit));

    %% ── Linear forced-origin fit ─────────────────────────────────────────
    C_lin    = Az_fit \ F_fit;
    pred_lin = C_lin .* Az_fit;
    r2_lin   = 1 - sum((F_fit-pred_lin).^2) / sum((F_fit-mean(F_fit)).^2);
    rmse_lin = sqrt(mean((F_fit-pred_lin).^2));
    K_eff    = C_lin * m_foot / (rho_s * g_cm);

    %% ── Power law fit (grid search p=0.20-1.50) ──────────────────────────
    p_grid = 0.20:0.01:1.50;
    best_p = struct('R2',-Inf,'p',NaN,'C',NaN,'rmse',Inf);
    for pp = p_grid
        Az_pow  = Az_fit.^pp;
        C_pp    = Az_pow \ F_fit;
        pred_pp = C_pp .* Az_pow;
        r2_pp   = 1 - sum((F_fit-pred_pp).^2)/sum((F_fit-mean(F_fit)).^2);
        if r2_pp > best_p.R2
            best_p.R2   = r2_pp;
            best_p.p    = pp;
            best_p.C    = C_pp;
            best_p.rmse = sqrt(mean((F_fit-pred_pp).^2));
        end
    end

    fprintf('\n-- plot_fz_vs_Az: LFFM geometry correction ------------------\n');
    fprintf('Power law: C=%.1f  p=%.2f  R2=%.4f  RMSE=%.1f\n', ...
        best_p.C, best_p.p, best_p.R2, best_p.rmse);
    fprintf('\nGeometry comparison:\n');
    fprintf('  Katsuragi sphere: k/m=%.0f  K_eff=%.2f  slope_Az=%.1f\n', ...
        km_sphere, K_eff_sphere, slope_B);
    fprintf('  Kang theory:      K_phi=%.1f  slope_Az=%.1f\n', ...
        K_phi_theory, slope_A);
    fprintf('  Jerboa foot:      k/m=%.0f  K_eff=%.1f  C_lin=%.1f\n', ...
        k_over_m, K_eff, C_lin);
    fprintf('  Foot/sphere K_eff ratio: %.1fx\n', K_eff/K_eff_sphere);
    fprintf('  Foot/theory K_eff ratio: %.1fx\n', K_eff/K_phi_theory);
    fprintf('-------------------------------------------------------------\n\n');

    %% ── Axis limits ──────────────────────────────────────────────────────
    Az_max = max(Az_fit) * 1.18;
    F_max  = max(F_fit + SE_fit) * 1.30;
    if ~isfinite(F_max) || F_max <= 0, F_max = max(F_fit)*1.30; end

    %% ── Figure ───────────────────────────────────────────────────────────
    fig = figure('Name','Step 3: F(z)/m vs A(z)*z', ...
                 'ToolBar','none','MenuBar','none');
    fig.Position = [220 80 800 620];
    ax = axes(fig, 'Position',[0.11 0.12 0.85 0.83]);
    hold(ax,'on');

    Az_line = linspace(0, Az_max, 300);

    % y=0 reference
    yline(ax, 0, '--', 'Color',[0.10 0.10 0.10], 'LineWidth',2.4, ...
        'HandleVisibility','off');

    % Kang theory — grey dotted (Approach A, quasistatic compact)
    plot(ax, Az_line, slope_A .* Az_line, ':', ...
        'Color',col_kang, 'LineWidth',4.0, ...
        'DisplayName', sprintf('Kang theory (K_phi=%.1f, compact)', K_phi_theory));

    % Katsuragi sphere — black dotted (Approach B, dynamic)
    plot(ax, Az_line, slope_B .* Az_line, ':', ...
        'Color',col_kats, 'LineWidth',4.0, ...
        'DisplayName', sprintf('Katsuragi sphere (k/m=%.0f s^{-2})', km_sphere));

    % Your foot — dark green dashed (power law)
    plot(ax, Az_line, best_p.C .* Az_line.^best_p.p, '--', ...
        'Color',col_green, 'LineWidth',3.5, ...
        'DisplayName', sprintf('F/m = %.0f [A(z)z]^{%.2f},  R^2=%.2f', ...
            best_p.C, best_p.p, best_p.R2));

    % Stable symbols — crimson-pink fill, black edge
    for ii = 1:numel(Az_fit)
        errorbar(ax, Az_fit(ii), F_fit(ii), SE_fit(ii), SE_fit(ii), 0, 0, ...
            'o', 'Color','k', 'MarkerFaceColor',col_stable, ...
            'MarkerEdgeColor','k', 'MarkerSize',ms, ...
            'LineWidth',lw_sym, 'CapSize',3, 'HandleVisibility','off');
    end


    set(ax,'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in', ...
        'XLim',[0, Az_max],'YLim',[0, F_max]);
    grid(ax,'off');
    xlabel(ax,'$A_{\rm hull}(z)\cdot z$  (cm$^3$)','FontSize',22, ...
        'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax,'$F(z)/m$  (cm s$^{-2}$)','FontSize',22, ...
        'Interpreter','latex','Color',[0 0 0]);


    %% ── Pack fitStats ────────────────────────────────────────────────────
    fitStats.C_lin         = C_lin;
    fitStats.r2_lin        = r2_lin;
    fitStats.rmse_lin      = rmse_lin;
    fitStats.K_eff         = K_eff;
    fitStats.C_pow         = best_p.C;
    fitStats.p_pow         = best_p.p;
    fitStats.r2_pow        = best_p.R2;
    fitStats.rmse_pow      = best_p.rmse;
    fitStats.slope_B       = slope_B;
    fitStats.slope_A       = slope_A;
    fitStats.K_eff_sphere  = K_eff_sphere;
    fitStats.K_phi_theory  = K_phi_theory;
    fitStats.km_sphere     = km_sphere;
    fitStats.A_hull        = A_hull;
    fitStats.Az            = Az;
    fitStats.Az_fit        = Az_fit;
    fitStats.z_min_fit     = z_min_fit;
    fitStats.d1            = d1;
    fitStats.k_over_m      = k_over_m;
end