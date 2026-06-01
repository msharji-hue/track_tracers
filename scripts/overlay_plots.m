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

%% ── FIGURE 4: Normalized collapse ────────────────────────────────────────
fig4 = plot_normalized_collapse(heights, cmap);

%% ── FIGURE 5: a+g vs v² — extracts d1_kinematic (Katsuragi Fig. 3a) ─────
z_targets = [1.0, 1.3, 1.6, 1.9, 2.2];

[fig5, d1_kinematic, Fz_over_m, fitStats] = ...
    plot_ag_vs_v2(heights, cmap, z_targets, 40);

fprintf('d1 = %.3f cm  |  R2 = %.4f\n', ...
    d1_kinematic, fitStats.r2);

%% ── FIGURE 6: d vs v0 — extracts d0 (Katsuragi Fig. 2b) ─────────────────
[fig6, d0] = plot_depth_vs_v0(heights, cmap);

fprintf('d0 = %.3f cm\n', d0);


%% ── FIGURE 7: F(z)/m vs z — extracts k_over_m (Katsuragi Fig. 3b) ───────
[fig7, k_over_m, fitStatsFz] = plot_fz_vs_z(heights, d1_kinematic, cmap, ...
    Fz_over_m, fitStats.z_targets, d0, fitStats.Fz_over_m_se);

fprintf('k/m = %.1f s^-2\n', k_over_m);
 
%% ── FIGURE 8: t_stop vs v0 — characteristic time scale ──────────────────
[fig8, stats8] = plot_tstop_vs_v0(heights, cmap, ...
    'stlFile', 'jerboa_foot_model_rectangularbeam.stl', ...
    'alpha', 1.5);

fprintf('A_ref = %.4f cm^2\n', stats8.A_ref_cm2);
fprintf('D_eff = %.3f cm\n', stats8.D_eff_cm);
fprintf('Vc    = %.2f cm/s\n', stats8.Vc_cm_s);
fprintf('Tc    = %.4f s\n', stats8.Tc_s);
 
%% ── CLEAN MEAN-LINE PLOTS ────────────────────────────────────────────────
fig_z_clean  = plot_mean_lines(heights, 'z_smooth',  '$z$  (cm)',             cmap);
fig_v_clean  = plot_mean_lines(heights, 'v_smooth',  '$v$  (cm s$^{-1}$)',   cmap, ...
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
    allFigs  = {fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, ...
                fig_z_clean, fig_v_clean, fig_ag_clean, fig_ex_z, fig_ex_v, fig_ex_ag};
    allNames = {'z_vs_t', 'v_vs_t', 'aplusg_vs_t', 'normalized_collapse', ...
                'ag_vs_v2', 'fz_vs_z', 'd_vs_v0', 'tstop_vs_v0', ...
                'logag_vs_logv', 'v2_vs_z', ...
                'z_vs_t_clean', 'v_vs_t_clean', 'aplusg_vs_t_clean','example_z', 'example_v', 'example_ag'};
    save_all_figures(allFigs, allNames, outDir);
    save_analysis(heights, hVals, cmap, outDir);
end