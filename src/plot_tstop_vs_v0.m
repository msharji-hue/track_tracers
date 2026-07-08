function [fig, stats] = plot_tstop_vs_v0(heights, cmap, varargin)
% PLOT_TSTOP_VS_V0  Stopping time vs impact speed.
%   Replicates Katsuragi & Durian 2007 Fig. 1b inset.
%
%   Key additions vs original:
%     1. Forward model t_stop prediction (black dotted) — requires d1, k_over_m
%     2. D_eff computed from stable-shaft A_hull only (z >= 1.30 cm)
%        matching the regime where d1 and k/m are calibrated
%
%   Characteristic scales:
%     Lc  = D_eff = sqrt(4*A_hull_stable_mean/pi)
%     Vc  = sqrt(Lc*g)  — velocity scale
%     Tc  = sqrt(Lc/g)  — time scale
%
%   For Katsuragi sphere: Lc = Db = 2.54 cm, Vc = 49.9 cm/s, Tc = 0.051 s
%
%   Usage:
%     [fig8, stats8] = plot_tstop_vs_v0(heights, cmap, ...
%         'footArea', footArea, 'd1', d1_kinematic, 'k_over_m', k_over_m);
%
%   Name-value options:
%     'footArea'       - struct from extract_foot_area_vs_depth
%     'stlFile'        - STL filename (fallback if footArea not passed)
%     'alpha'          - alpha for STL extraction
%     'd1'             - inertial length scale [cm] for forward model
%     'k_over_m'       - friction coefficient [s^-2] for forward model
%     'showScaleLines' - true/false (default true)
%     'useSeconds'     - true/false (default true)
%     'z_min_stable'   - stable shaft boundary [cm] (default 1.30)

    %% ── Parse inputs ─────────────────────────────────────────────────────
    p = inputParser;
    addParameter(p,'footArea',      [],     @(x) isempty(x)||isstruct(x));
    addParameter(p,'stlFile', ...
        'jerboa_foot_model_rectangularbeam.stl', @(s) ischar(s)||isstring(s));
    addParameter(p,'alpha',         1.5,    @isnumeric);
    addParameter(p,'d1',            [],     @isnumeric);
    addParameter(p,'k_over_m',      [],     @isnumeric);
    addParameter(p,'showScaleLines',true,   @islogical);
    addParameter(p,'useSeconds',    true,   @islogical);
    addParameter(p,'z_min_stable',  1.30,   @isnumeric);
    parse(p, varargin{:});
    opt = p.Results;

    g_cm         = 980;
    nH           = numel(heights);
    doForward    = ~isempty(opt.d1) && ~isempty(opt.k_over_m);
    z_min_stable = opt.z_min_stable;

    %% ── Get foot area profile ────────────────────────────────────────────
    footArea = opt.footArea;
    if isempty(footArea)
        fprintf('\nExtracting STL area profile...\n');
        footArea = extract_foot_area_vs_depth(char(opt.stlFile), ...
            'alpha', opt.alpha);
    end
    if ~isfield(footArea,'depth_cm') || ~isfield(footArea,'A_hull_sm')
        error('footArea must contain depth_cm and A_hull_sm.');
    end

    %% ── D_eff from stable-shaft A_hull only ──────────────────────────────
    % Use z >= z_min_stable — same regime where d1 and k/m are calibrated
    in_stable = footArea.depth_cm >= z_min_stable & ...
                isfinite(footArea.A_hull_sm) & ...
                footArea.A_hull_sm > 0;

    if ~any(in_stable)
        warning('No stable-shaft depths found in footArea. Using full range.');
        in_stable = isfinite(footArea.A_hull_sm) & footArea.A_hull_sm > 0;
    end

    A_ref_cm2 = mean(footArea.A_hull_sm(in_stable), 'omitnan');
    D_eff_cm  = sqrt(4*A_ref_cm2/pi);
    Vc_cm_s   = sqrt(D_eff_cm * g_cm);
    Tc_s      = sqrt(D_eff_cm / g_cm);

    fprintf('\n-- t_stop vs v0 ---------------------------------------------\n');
    fprintf('D_eff (stable shaft) = %.4f cm  (from z >= %.2f cm)\n', ...
        D_eff_cm, z_min_stable);
    fprintf('Vc = sqrt(D_eff*g)   = %.2f cm/s\n', Vc_cm_s);
    fprintf('Tc = sqrt(D_eff/g)   = %.5f s\n', Tc_s);

    % Compare to Katsuragi sphere
    D_sphere = 2.54;
    fprintf('Katsuragi sphere: Db=%.2f cm, Vc=%.1f cm/s, Tc=%.4f s\n', ...
        D_sphere, sqrt(D_sphere*g_cm), sqrt(D_sphere/g_cm));
    fprintf('D_eff/Db = %.3f\n', D_eff_cm/D_sphere);
    fprintf('-------------------------------------------------------------\n\n');

    %% ── Collect plot data ────────────────────────────────────────────────
    v0_all     = [];
    tstop_all  = [];
    for j = 1:nH
        v0_all    = [v0_all,    heights(j).v0_mean];
        tstop_all = [tstop_all, heights(j).tstop_mean];
    end

    if opt.useSeconds
        t_scale = 1;
        y_label = '$t_{\rm stop}$  (s)';
        Tc_plot = Tc_s;
    else
        t_scale = 1000;
        y_label = '$t_{\rm stop}$  (ms)';
        Tc_plot = Tc_s * 1000;
    end

    %% ── Forward model t_stop prediction ─────────────────────────────────
    v0_line   = linspace(1, max(v0_all)*1.08, 150);
    t_fwd     = nan(size(v0_line));

    if doForward
        d1_val   = opt.d1;
        km_val   = opt.k_over_m;
        dt       = 0.00005;   % s — fine timestep for accuracy
        max_time = 0.5;       % s — safety cap

        for vi = 1:numel(v0_line)
            v_t = v0_line(vi);
            z_t = 0;
            t_t = 0;
            while v_t > 0 && t_t < max_time
                d1_z  = max(d1_val, 3.40 - 1.01*z_t);   % depth-varying d1
                accel = g_cm - km_val*z_t - v_t^2/d1_z;
                v_new = v_t + accel*dt;
                if v_new <= 0
                    % Interpolate stopping time
                    t_t = t_t + dt * v_t / (v_t - v_new);
                    break;
                end
                v_t = v_new;
                z_t = z_t + v_t*dt;
                t_t = t_t + dt;
            end
            t_fwd(vi) = t_t;
        end

        if opt.useSeconds
            t_fwd_plot = t_fwd;
        else
            t_fwd_plot = t_fwd * 1000;
        end
    end

    %% ── Figure ───────────────────────────────────────────────────────────
    fig = figure('Name','t_stop vs v0','ToolBar','none','MenuBar','none');
    fig.Position = [100 100 580 480];
    ax = axes(fig,'Position',[0.13 0.13 0.83 0.82]);
    hold(ax,'on');

    % Forward model — black dotted (plot first so data is on top)
    if doForward
        plot(ax, v0_line, t_fwd_plot, 'k:', 'LineWidth',2.2, ...
            'HandleVisibility','off');
    end

    % Data points per height group
    for j = 1:nH
        hg = heights(j);
        [~, marker] = get_height_style(hg.h_cm);
        errorbar(ax, hg.v0_mean, hg.tstop_mean*t_scale, ...
            hg.tstop_std*t_scale, hg.tstop_std*t_scale, ...
            hg.v0_std, hg.v0_std, marker, ...
            'Color',cmap(j,:), 'MarkerFaceColor','none', ...
            'MarkerEdgeColor',cmap(j,:), 'MarkerSize',9, 'LineWidth',1.8, ...
            'HandleVisibility','off');
    end

    % Characteristic scale guide lines
    if opt.showScaleLines
        xline(ax, Vc_cm_s, '--', 'Color',[0.45 0.45 0.45], ...
            'LineWidth',1.4,'HandleVisibility','off');
        yline(ax, Tc_plot,  '--', 'Color',[0.45 0.45 0.45], ...
            'LineWidth',1.4,'HandleVisibility','off');

        % Labels on guide lines
        text(ax, Vc_cm_s + 2, max(tstop_all*t_scale)*0.15, ...
            sprintf('$V_c=%.1f$ cm s$^{-1}$', Vc_cm_s), ...
            'FontSize',9,'Interpreter','latex', ...
            'Color',[0.35 0.35 0.35], ...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
        text(ax, max(v0_all)*1.06, Tc_plot, ...
            sprintf('$T_c=%.4f$ s', Tc_s), ...
            'FontSize',9,'Interpreter','latex', ...
            'Color',[0.35 0.35 0.35], ...
            'HorizontalAlignment','right','VerticalAlignment','bottom');
    end

    set(ax,'FontSize',13,'Box','on','LineWidth',1.2, ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XMinorTick','on','YMinorTick','on','TickDir','in', ...
        'XLim',[0, max(v0_all)*1.10], ...
        'YLim',[0, max(tstop_all*t_scale)*1.28]);
    grid(ax,'off');
    xlabel(ax,'$v_0$  (cm s$^{-1}$)','FontSize',16, ...
        'Interpreter','latex','Color',[0 0 0]);
    ylabel(ax, y_label,'FontSize',16, ...
        'Interpreter','latex','Color',[0 0 0]);

    %% ── Annotation ───────────────────────────────────────────────────────
    ann_lines = {sprintf('$L_c = D_{\\rm eff} = %.3f$ cm', D_eff_cm), ...
                 sprintf('$V_c = \\sqrt{L_c g} = %.1f$ cm s$^{-1}$', Vc_cm_s), ...
                 sprintf('$T_c = \\sqrt{L_c/g} = %.4f$ s', Tc_s), ...
                 '(stable shaft, $z \geq 1.30$ cm)'};
    if doForward
        ann_lines{end+1} = sprintf(...
            'Fwd model: $d_1=%.2f$ cm, $k/m=%.0f$ s$^{-2}$', ...
            opt.d1, opt.k_over_m);
    end

    text(ax, 0.97, 0.97, strjoin(ann_lines, '\n'), ...
        'Units','normalized','Interpreter','latex','FontSize',9, ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25], ...
        'Margin',4);

    %% ── Output stats ─────────────────────────────────────────────────────
    stats.A_ref_cm2   = A_ref_cm2;
    stats.D_eff_cm    = D_eff_cm;
    stats.Lc_cm       = D_eff_cm;
    stats.Vc_cm_s     = Vc_cm_s;
    stats.Tc_s        = Tc_s;
    stats.g_cm_s2     = g_cm;
    stats.z_min_stable= z_min_stable;
    if doForward
        stats.t_fwd   = t_fwd;
        stats.v0_line = v0_line;
    end
end