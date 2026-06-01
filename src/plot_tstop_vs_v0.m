function [fig, stats] = plot_tstop_vs_v0(heights, cmap, varargin)
% PLOT_TSTOP_VS_V0  Stopping time vs impact speed.
%
%   Matches the style/logic of Katsuragi & Durian Fig. 1b inset:
%       t_stop vs v0
%       characteristic velocity scale sqrt(D_eff*g)
%       characteristic time scale sqrt(D_eff/g)
%
%   For a sphere, Katsuragi uses projectile diameter Db.
%   For the jerboa-foot model, this function estimates an equivalent
%   characteristic size from the bare printed STL area:
%
%       D_eff = sqrt(4*A_bare_ref/pi)
%
%   Usage:
%       [fig8, stats8] = plot_tstop_vs_v0(heights, cmap);
%
%       [fig8, stats8] = plot_tstop_vs_v0(heights, cmap, ...
%           'stlFile', 'jerboa_foot_model_rectangularbeam.stl');
%
%       [fig8, stats8] = plot_tstop_vs_v0(heights, cmap, ...
%           'footArea', out);
%
%   Name-value options:
%       'footArea'       - output from extract_foot_area_vs_depth
%       'stlFile'        - STL filename used if footArea is not passed
%       'alpha'          - alpha value for STL extraction fallback
%       'showScaleLines' - true/false for Vc and Tc guide lines
%       'useSeconds'     - true/false, true matches Katsuragi inset
%
%   Outputs:
%       fig   - figure handle
%       stats - struct with A_ref, D_eff, Vc, Tc, etc.

    %% ── Parse inputs ─────────────────────────────────────────────────────
    p = inputParser;

    addParameter(p, 'footArea', [], @(x) isempty(x) || isstruct(x));
    addParameter(p, 'stlFile', 'jerboa_foot_model_rectangularbeam.stl', ...
        @(s) ischar(s) || isstring(s));
    addParameter(p, 'alpha', 1.5, @isnumeric);
    addParameter(p, 'showScaleLines', true, @islogical);
    addParameter(p, 'useSeconds', true, @islogical);

    parse(p, varargin{:});
    opt = p.Results;

    g_cm = 980;
    nH   = numel(heights);

    %% ── Get or extract bare STL area profile ─────────────────────────────
    footArea = opt.footArea;

    if isempty(footArea)
        fprintf('\nExtracting STL area profile for D_eff calculation...\n');

        footArea = extract_foot_area_vs_depth(char(opt.stlFile), ...
            'alpha', opt.alpha);

        fprintf('Finished STL area extraction for stopping-time scale.\n\n');
    end

    if ~isfield(footArea, 'depth_cm') || ~isfield(footArea, 'A_bare_sm')
        error('footArea must contain fields depth_cm and A_bare_sm.');
    end

    %% ── Determine experimental penetration range ─────────────────────────
    max_d_exp = 0;

    for j = 1:nH
        if isfield(heights(j), 'd_mean')
            max_d_exp = max(max_d_exp, max(heights(j).d_mean, [], 'omitnan'));
        else
            for i = 1:heights(j).nTrials
                max_d_exp = max(max_d_exp, ...
                    heights(j).trials(i).scalars.d_final_cm);
            end
        end
    end

    in_range = footArea.depth_cm >= 0 & ...
               footArea.depth_cm <= 1.10*max_d_exp & ...
               isfinite(footArea.A_bare_sm);

    if ~any(in_range)
        error('No valid STL area values found within the experimental penetration range.');
    end

%% ── Characteristic size from convex-hull intrusion envelope ───────────
% Katsuragi uses sphere diameter Db.
% For the foot, use an equivalent diameter from the convex hull area,
% which better represents the swept envelope seen by the grains.

if isfield(footArea, 'A_hull_sm')
    A_ref_cm2 = mean(footArea.A_hull_sm(in_range), 'omitnan');
elseif isfield(footArea, 'A_hull')
    A_ref_cm2 = mean(footArea.A_hull(in_range), 'omitnan');
else
    error('footArea must contain A_hull_sm or A_hull for convex-hull scaling.');
end

D_eff_cm = sqrt(4*A_ref_cm2/pi);

% Optional characteristic length override
% Leave empty to use convex-hull equivalent diameter.
% Try values like 2.0, 4.0, or 6.5 cm to see how both guide lines move.
Lc_override_cm = [];

if isempty(Lc_override_cm)
    Lc_cm = D_eff_cm;
    Lc_label = 'convex-hull envelope';
else
    Lc_cm = Lc_override_cm;
    Lc_label = 'manual length scale';
end

Vc_cm_s = sqrt(Lc_cm*g_cm);
Tc_s    = sqrt(Lc_cm/g_cm);

    %% ── Collect plot data ────────────────────────────────────────────────
    v0_mean_all    = [];
    tstop_mean_all = [];

    for j = 1:nH
        v0_mean_all    = [v0_mean_all, heights(j).v0_mean];
        tstop_mean_all = [tstop_mean_all, heights(j).tstop_mean];
    end

    if opt.useSeconds
        t_scale = 1;
        y_label = '$t_\mathrm{stop}$  (s)';
        Tc_plot = Tc_s;
    else
        t_scale = 1000;
        y_label = '$t_\mathrm{stop}$  (ms)';
        Tc_plot = Tc_s * 1000;
    end

    %% ── Figure ───────────────────────────────────────────────────────────
    fig = figure('Name','t_stop vs v0', ...
                 'ToolBar','none', ...
                 'MenuBar','none');

    fig.Position = [100 100 540 440];

    ax = axes(fig, 'Position', [0.14 0.13 0.78 0.82]);
    hold(ax,'on');

    for j = 1:nH

        hg = heights(j);

        [~, marker] = get_height_style(hg.h_cm);

        errorbar(ax, hg.v0_mean, hg.tstop_mean*t_scale, ...
            hg.tstop_std*t_scale, hg.tstop_std*t_scale, ...
            hg.v0_std, hg.v0_std, marker, ...
            'Color', cmap(j,:), ...
            'MarkerFaceColor','none', ...
            'MarkerEdgeColor', cmap(j,:), ...
            'MarkerSize', 9, ...
            'LineWidth', 1.8, ...
            'HandleVisibility','off');
    end

    %% ── Characteristic scale guide lines ─────────────────────────────────
    if opt.showScaleLines

        xline(ax, Vc_cm_s, '--', ...
            'Color', [0.45 0.45 0.45], ...
            'LineWidth', 1.4, ...
            'HandleVisibility','off');

        yline(ax, Tc_plot, '--', ...
            'Color', [0.45 0.45 0.45], ...
            'LineWidth', 1.4, ...
            'HandleVisibility','off');
    end

    %% ── Axes styling ─────────────────────────────────────────────────────
    set(ax,'FontSize',13, ...
        'Box','on', ...
        'LineWidth',1.2, ...
        'XColor',[0 0 0], ...
        'YColor',[0 0 0], ...
        'XMinorTick','on', ...
        'YMinorTick','on', ...
        'TickDir','in');

    grid(ax,'off');

    xlabel(ax,'$v_0$  (cm s$^{-1}$)', ...
        'FontSize',16, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);

    ylabel(ax, y_label, ...
        'FontSize',16, ...
        'Interpreter','latex', ...
        'Color',[0 0 0]);

    xlim(ax, [0, max(v0_mean_all)*1.08]);
    ylim(ax, [0, max(tstop_mean_all*t_scale)*1.25]);

    %% ── Annotation ───────────────────────────────────────────────────────
    if opt.useSeconds
        txt = sprintf(['$D_{\\rm eff}=%.2f$ cm\n' ...
                       '$\\sqrt{D_{\\rm eff}g}=%.1f$ cm s$^{-1}$\n' ...
                       '$\\sqrt{D_{\\rm eff}/g}=%.3f$ s'], ...
            D_eff_cm, Vc_cm_s, Tc_s);
    else
        txt = sprintf(['$D_{\\rm eff}=%.2f$ cm\n' ...
                       '$\\sqrt{D_{\\rm eff}g}=%.1f$ cm s$^{-1}$\n' ...
                       '$\\sqrt{D_{\\rm eff}/g}=%.1f$ ms'], ...
            D_eff_cm, Vc_cm_s, Tc_s*1000);
    end

    text(ax, 0.97, 0.95, ...
    sprintf(['$L_c=%.2f$ cm\n' ...
             '$V_c=\\sqrt{L_c g}=%.1f$ cm s$^{-1}$\n' ...
             '$T_c=\\sqrt{L_c/g}=%.3f$ s\n' ...
             '%s'], ...
        Lc_cm, Vc_cm_s, Tc_s, Lc_label), ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'FontSize',10.5, ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'BackgroundColor',[1 1 1], ...
    'EdgeColor',[0.25 0.25 0.25], ...
    'Margin',4);

    %% ── Output stats ─────────────────────────────────────────────────────
stats.A_ref_cm2 = A_ref_cm2;
stats.D_eff_cm  = D_eff_cm;
stats.Lc_cm     = Lc_cm;
stats.Vc_cm_s   = Vc_cm_s;
stats.Tc_s      = Tc_s;
stats.g_cm_s2   = g_cm;
end