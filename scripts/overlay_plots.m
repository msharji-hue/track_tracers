%% scripts/overlay_plots.m
clear; close all; clc;

codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

%% ── LOAD + GROUP ─────────────────────────────────────────────────────────
trials  = load_selected_trials();
%%
heights = group_trials_by_height(trials);
hVals   = [heights.h_cm];
nH      = numel(heights);
cmap    = turbo(nH);

%% ── OPTIONAL REMOVAL ─────────────────────────────────────────────────────
fprintf('\n── Height groups ─────────────────────────────────\n');
for j = 1:nH
    fprintf('  [%d]  h=%5.2f cm | n=%d | v0=%5.1f±%4.1f cm/s | d=%.3f±%.3f cm\n', ...
        j, heights(j).h_cm, heights(j).nTrials, ...
        heights(j).v0_mean, heights(j).v0_std, ...
        heights(j).d_mean,  heights(j).d_std);
end
removeGroups = input('Remove height groups by index (e.g. [2 5]), or []: ');
if ~isempty(removeGroups)
    heights(removeGroups) = [];
    hVals = [heights.h_cm];
    nH    = numel(heights);
    cmap  = turbo(nH);
    fprintf('Removed %d group(s). %d remaining.\n', numel(removeGroups), nH);
end

%% ── FIGURE 1: z vs t ─────────────────────────────────────────────────────
[ax1, ~] = plot_time_series(heights, 'z_smooth', 'd_final_cm', 1, ...
    'z vs t', '$z$  (cm)', 'none', cmap);
legend(ax1, 'off');
yline(ax1, 0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
fig1 = ax1.Parent;

%% ── FIGURE 2: v vs t ─────────────────────────────────────────────────────
[ax2, ~] = plot_time_series(heights, 'v_smooth', 'v0_cm_s', 1, ...
    'v vs t', '$v$  (cm s$^{-1}$)', 'none', cmap);
legend(ax2, 'off');
yline(ax2, 0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'HandleVisibility', 'off');
for j = 1:nH
    [~, marker] = get_height_style(heights(j).h_cm);
    errorbar(ax2, 0, heights(j).v0_mean, heights(j).v0_std, marker, ...
        'Color', cmap(j,:), 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', cmap(j,:), 'MarkerSize', 9, ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');
end
fig2 = ax2.Parent;
ylim(ax2, [0, ax2.YLim(2)]);

%% ── FIGURE 3: a+g vs t ───────────────────────────────────────────────────
fig3 = figure('Name','a+g vs t','ToolBar','none','MenuBar','none');
fig3.Position = [100 100 980 520];
ax3 = axes(fig3, 'Position', [0.10 0.13 0.74 0.82]);
hold(ax3, 'on');
for j = 1:nH
    for i = 1:heights(j).nTrials
        k        = heights(j).trials(i).kinematics;
        idx_plot = k.impact_index:5:k.stopFrame;
        plot(ax3, k.t_s(idx_plot), k.a_plus_g(idx_plot), '.', ...
            'Color', [cmap(j,:), 0.20], 'MarkerSize', 4, ...
            'LineStyle', 'none', 'HandleVisibility', 'off');
    end
end
for j = 1:nH
    hg = heights(j);
    [tPhys, yMean, yStd, nValid] = interp_to_norm_time(hg, 'a_plus_g');
    lw = 2.5; if hg.nTrials < 3, lw = 1.2; end
    plot_shaded_band(ax3, tPhys, yMean, yStd, nValid, cmap(j,:), '', 1, lw);
end
set(ax3, 'FontSize',13,'Box','on','LineWidth',1.2, ...
    'XColor',[0 0 0],'YColor',[0 0 0], ...
    'XMinorTick','on','YMinorTick','on','TickDir','in', ...
    'XLim',[0 ax3.XLim(2)]);
yline(ax3, 0, 'k-', 'LineWidth', 1.8, 'HandleVisibility', 'off');
ylim([0 ax3.YLim(2)]);
grid(ax3,'off');
xlabel(ax3,'$t$  (s)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
ylabel(ax3,'$a+g$  (cm s$^{-2}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);

%% ── FIGURE 5 + Steps 1/2/3 shared computation ───────────────────────────
z_targets_all = (0.10 : 0.10 : 2.2)';
z_display     = [1.3, 1.5, 1.7, 1.9, 2.1]';

%% ── FIGURE 4a: Normalized collapse (linear axes) ────────────────────────
fig4a = plot_normalized_collapse(heights, cmap);

%% ── FIGURE 4b: Power-law proof (log-log) ────────────────────────────────
% Runs after d1 and k_over_m are available — place after fig7 call.
% Declared here; called after fig7 below.

%% ── FIGURE 5: a+g vs v ───────────────────────────────────────────────────
[fig5, d1_kinematic, Fz_over_m, fitStats5] = ...
    plot_ag_vs_v2(heights, cmap, z_targets_all, z_display, 40);

fprintf('d1 = %.3f cm  (CV = %.1f%%)\n', d1_kinematic, fitStats5.d1_cv);

%% ── FIGURE 6: d vs v0 ────────────────────────────────────────────────────
[fig6, d0] = plot_depth_vs_v0(heights, cmap, [], d1_kinematic);
fprintf('d0 = %.3f cm\n', d0);

%% ── Load STL foot area ───────────────────────────────────────────────────
stlMatFile = fullfile('/Users/muhannadalsharji/Documents/track_tracers/data', ...
    'jerboa_foot_model_rectangularbeam_area_profile.mat');
if exist(stlMatFile, 'file')
    tmp      = load(stlMatFile, 'out');
    footArea = tmp.out;
    fprintf('Loaded STL area profile: %d depth levels\n', numel(footArea.depth_cm));
else
    fprintf('STL area profile not found — Steps 2 and 3 figures skipped.\n');
    footArea = [];
end

%% ── FIGURE 7  (Step 1): F(z)/m vs z ─────────────────────────────────────
[fig7, k_over_m, fitStats7] = plot_fz_vs_z( ...
    heights, d1_kinematic, cmap, ...
    z_targets_all, fitStats5.Fz_over_m, ...
    fitStats5.Fz_per_height, fitStats5.Fz_per_height_se, ...
    fitStats5.v0_per_height, d0);

fprintf('k/m = %.1f s^-2  |  R2 = %.4f\n', k_over_m, fitStats7.r2_lin);

%% ── FIGURE 4b: Power-law collapse — placed here so d1 and k/m are ready ──
fig4b = plot_power_law_collapse(heights, cmap, d1_kinematic, k_over_m);

%% ── FIGURE 7b (Step 2): F(z)/m vs V_hull(z) ─────────────────────────────
% SE_plot from fitStats7 provides OLS-based SE at all 22 depths —
% ensures pre-stable symbols have vertical error bars too.
fig7b = [];
if ~isempty(footArea)
    [fig7b, fitStats7b] = plot_fz_vs_V( ...
        z_targets_all, fitStats5.Fz_over_m, fitStats7.SE_plot, ...
        d1_kinematic, footArea);
    fprintf('Kang slope = %.1f  K_eff = %.1f  R2 = %.4f\n', ...
        fitStats7b.slope_fit, fitStats7b.K_eff, fitStats7b.r2_fit);
end

%% ── FIGURE 7c (Step 3): F(z)/m vs A(z)*z ────────────────────────────────
fig7c = [];
if ~isempty(footArea)
    [fig7c, fitStats7c] = plot_fz_vs_Az( ...
        z_targets_all, fitStats5.Fz_over_m, fitStats7.SE_plot, ...
        d1_kinematic, k_over_m, footArea);
    fprintf('A(z)*z linear: C=%.1f  K_eff=%.1f  R2=%.4f\n', ...
        fitStats7c.C_lin, fitStats7c.K_eff, fitStats7c.r2_lin);
    fprintf('A(z)*z power:  C=%.1f  p=%.2f  R2=%.4f\n', ...
        fitStats7c.C_pow, fitStats7c.p_pow, fitStats7c.r2_pow);
end

%% ── FIGURE 6b: d vs v0 with forward model overlay ───────────────────────
[fig6b, ~] = plot_depth_vs_v0(heights, cmap, [], d1_kinematic, k_over_m);
%% ── FIGURE 8: t_stop vs v0 — characteristic time scale ──────────────────
[fig8, stats8] = plot_tstop_vs_v0(heights, cmap, ...
    'stlFile', 'jerboa_foot_model_rectangularbeam.stl', ...
    'alpha', 1.5);

fprintf('A_ref = %.4f cm^2\n', stats8.A_ref_cm2);
fprintf('D_eff = %.3f cm  (convex-hull equivalent diameter)\n', stats8.D_eff_cm);
fprintf('Lc    = %.3f cm  (length scale used for guide lines)\n', stats8.Lc_cm);
fprintf('Vc    = %.2f cm/s\n', stats8.Vc_cm_s);
fprintf('Tc    = %.4f s\n', stats8.Tc_s);

%% ── CLEAN MEAN-LINE PLOTS ────────────────────────────────────────────────
fig_z_clean  = plot_mean_lines(heights, 'z_smooth',  '$z$  (cm)',           cmap);
fig_v_clean  = plot_mean_lines(heights, 'v_smooth',  '$v$  (cm s$^{-1}$)', cmap, ...
    'yline', true, 'v0marks', true);
fig_ag_clean = plot_mean_lines(heights, 'a_plus_g',  '$a+g$  (cm s$^{-2}$)', cmap, ...
    'xlim0', true, 'ylim0', true);
for f = [fig_z_clean fig_v_clean fig_ag_clean]
    legend(f.CurrentAxes, 'off');
end

[fig_ex_z, fig_ex_v, fig_ex_ag] = plot_example_trials(heights, cmap);

%% ── SAVE ─────────────────────────────────────────────────────────────────
outDir = uigetdir(pwd, 'Select folder to save figures');
if outDir ~= 0
    % Collect valid figure handles (fig7b may be empty if footArea missing)
    allFigs  = {fig1, fig2, fig3, fig4, fig5, fig6, fig6b, fig7};
    allNames = {'z_vs_t', 'v_vs_t', 'aplusg_vs_t', 'normalized_collapse', ...
                'ag_vs_v2', 'd_vs_v0', 'd_vs_v0_fwd', 'fz_vs_z'};

    if ~isempty(fig7b)
        allFigs{end+1}  = fig7b;
        allNames{end+1} = 'fz_vs_V';
    end

    allFigs  = [allFigs,  {fig_z_clean, fig_v_clean, fig_ag_clean, ...
                            fig_ex_z, fig_ex_v, fig_ex_ag}];
    allNames = [allNames, {'z_vs_t_clean', 'v_vs_t_clean', 'aplusg_vs_t_clean', ...
                            'example_z', 'example_v', 'example_ag'}];

    save_all_figures(allFigs, allNames, outDir);
    save_analysis(heights, hVals, cmap, outDir);
end