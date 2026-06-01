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

%% ── FIGURE 5: a+g vs v² — extracts d1_kinematic ─────────────────────────
z_targets = [1.0, 1.3, 1.6, 1.9, 2.2];
[fig5, d1_kinematic, ~, fitStats] = plot_ag_vs_v2(heights, cmap, z_targets, 40);

%% ── FIGURE 6: F(z)/m vs z ────────────────────────────────────────────────
[fig6, k_over_m] = plot_fz_vs_z(heights, d1_kinematic, cmap, ...
                                  fitStats.Fz_over_m, fitStats.z_targets);

fprintf('k/m = %.1f s^-2\n', k_over_m);
%% ── STL area extraction (needed for Fig 7) ───────────────────────────────
out = extract_foot_area_vs_depth('jerboa_foot_model_rectangularbeam.stl', ...
    'alpha', 1.5, 'gamma', 0.3);

%% ── FIGURE 7: d vs v0 with forward model ────────────────────────────────
fig7 = plot_depth_vs_v0(heights, cmap, out, d1_kinematic, k_over_m);

%% ── FIGURE 8: t_stop vs v0 ──────────────────────────────────────────────
fig8 = figure('Name','t_stop vs v0','ToolBar','none','MenuBar','none');
fig8.Position = [100 100 540 440];
ax8 = axes(fig8, 'Position', [0.14 0.13 0.78 0.82]);
hold(ax8, 'on');
for j = 1:nH
    hg = heights(j);
    [~, marker] = get_height_style(hg.h_cm);
    errorbar(ax8, hg.v0_mean, hg.tstop_mean*1000, ...
        hg.tstop_std*1000, hg.tstop_std*1000, hg.v0_std, hg.v0_std, marker, ...
        'Color', cmap(j,:), 'MarkerFaceColor','none', ...
        'MarkerEdgeColor', cmap(j,:), 'MarkerSize', 9, 'LineWidth', 1.8, ...
        'HandleVisibility', 'off');
end
set(ax8,'FontSize',13,'Box','on','LineWidth',1.2, ...
    'XColor',[0 0 0],'YColor',[0 0 0], ...
    'XMinorTick','on','YMinorTick','on','TickDir','in');
grid(ax8,'off');
xlabel(ax8,'$v_0$  (cm s$^{-1}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
ylabel(ax8,'$t_\mathrm{stop}$  (ms)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
add_h_colorbar(ax8, hVals, cmap);

%% ── FIGURE 9: log(a+g) vs log(v) ────────────────────────────────────────
fig9 = figure('Name','log(a+g) vs log(v)','ToolBar','none','MenuBar','none');
fig9.Position = [100 100 580 460];
ax9 = axes(fig9, 'Position', [0.14 0.13 0.76 0.82]);
hold(ax9, 'on');
for j = 1:nH
    for i = 1:heights(j).nTrials
        k   = heights(j).trials(i).kinematics;
        idx = k.impact_index:k.stopFrame;
        v   = k.v_smooth(idx);
        ag  = k.a_plus_g(idx);
        ok  = isfinite(v) & isfinite(ag) & v > 0 & ag > 0;
        plot(ax9, v(ok), ag(ok), '.', 'Color', [cmap(j,:), 0.40], ...
            'MarkerSize', 4, 'HandleVisibility', 'off');
    end
end
v_ref = logspace(0, log10(max([heights.v0_mean])*1.1), 100);
plot(ax9, v_ref, v_ref.^2 / d1_kinematic, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', '$v^2/d_1$  [inertial only]');
set(ax9,'XScale','log','YScale','log','FontSize',13,'Box','on', ...
    'LineWidth',1.2,'XColor',[0 0 0],'YColor',[0 0 0], ...
    'XMinorTick','on','YMinorTick','on','TickDir','in');
grid(ax9,'off');
xlabel(ax9,'$v$  (cm s$^{-1}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
ylabel(ax9,'$a+g$  (cm s$^{-2}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
legend(ax9,'show','FontSize',11,'Interpreter','latex','Box','off','Location','northwest');
text(ax9, 0.05, 0.92, sprintf('$d_1 = %.2f$ cm', d1_kinematic), ...
    'Units','normalized','Interpreter','latex','FontSize',11, ...
    'VerticalAlignment','top','Color',[0.15 0.15 0.15], ...
    'BackgroundColor',[1 1 1],'EdgeColor',[0.25 0.25 0.25],'Margin',4);
text(ax9, 0.97, 0.06, {'slope $= 2$ $\Rightarrow$ pure $v^2$ drag', ...
    'slope $< 2$ $\Rightarrow$ friction contributes'}, ...
    'Units','normalized','Interpreter','latex','FontSize',10, ...
    'HorizontalAlignment','right','Color',[0.25 0.25 0.25]);
add_h_colorbar(ax9, hVals, cmap);

%% ── FIGURE 10: v² vs z semi-log ─────────────────────────────────────────
fig10 = figure('Name','v^2 vs z','ToolBar','none','MenuBar','none');
fig10.Position = [100 100 580 460];
ax10 = axes(fig10, 'Position', [0.14 0.13 0.76 0.82]);
hold(ax10, 'on');
for j = 1:nH
    for i = 1:heights(j).nTrials
        k   = heights(j).trials(i).kinematics;
        idx = k.impact_index:k.stopFrame;
        z   = k.z_smooth(idx);
        v2  = k.v_smooth(idx).^2;
        ok  = isfinite(z) & isfinite(v2) & v2 > 0;
        plot(ax10, z(ok), v2(ok), '-', 'Color', [cmap(j,:), 0.55], 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end
end
set(ax10,'YScale','log','FontSize',13,'Box','on','LineWidth',1.2, ...
    'XColor',[0 0 0],'YColor',[0 0 0], ...
    'XMinorTick','on','YMinorTick','on','TickDir','in');
grid(ax10,'off');
xlabel(ax10,'$z$  (cm)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
ylabel(ax10,'$v^2$  (cm$^2$ s$^{-2}$)','FontSize',16,'Interpreter','latex','Color',[0 0 0]);
text(ax10, 0.05, 0.10, {'Straight line $=$ pure inertial decay', ...
    sprintf('slope $\\approx -2/d_1 = %.2f$ cm$^{-1}$', 2/d1_kinematic)}, ...
    'Units','normalized','Interpreter','latex','FontSize',10, ...
    'Color',[0.25 0.25 0.25],'VerticalAlignment','bottom');
add_h_colorbar(ax10, hVals, cmap);

%% ── CLEAN MEAN-LINE PLOTS ────────────────────────────────────────────────
fig_z_clean  = plot_mean_lines(heights, 'z_smooth',  '$z$  (cm)',             cmap);
fig_v_clean  = plot_mean_lines(heights, 'v_smooth',  '$v$  (cm s$^{-1}$)',   cmap, ...
                   'yline', true, 'v0marks', true);
fig_ag_clean = plot_mean_lines(heights, 'a_plus_g',  '$a+g$  (cm s$^{-2}$)', cmap, ...
                   'xlim0', true, 'ylim0', true);
for f = [fig_z_clean fig_v_clean fig_ag_clean]
    legend(f.CurrentAxes, 'off');
end

%% ── EXAMPLE TRIALS: low / mid / high v0 ─────────────────────────────────
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