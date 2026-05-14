%% scripts/overlay_plots.m
clear; close all; clc;

codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

%% LOAD + GROUP
trials  = load_selected_trials();
heights = group_trials_by_height(trials);
hVals   = [heights.h_cm];
nH      = numel(heights);
cmap    = parula(nH);

%% ── FIGURE 1: z vs t ─────────────────────────────────────────────────────
[ax1, ~] = plot_time_series(heights, 'z_smooth', 'd_final_cm', 10, ...
    'z vs t', '$z$  (cm)', 'none', cmap);
lgd1 = legend(ax1, 'show', 'FontSize', 10, 'Box', 'off');
lgd1.Location = 'none';
lgd1.Position = [0.87 0.30 0.12 0.40];

%% ── FIGURE 2: v vs t ─────────────────────────────────────────────────────
[ax2, ~] = plot_time_series(heights, 'v_smooth', 'v0_cm_s', 1, ...
    'v vs t', '$v$  (cm s$^{-1}$)', 'none', cmap);
yline(ax2, 0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'HandleVisibility', 'off');
for j = 1:nH
    hg          = heights(j);
    [~, marker] = get_height_style(hg.h_cm);
    errorbar(ax2, 0, hg.v0_mean, hg.v0_std, marker, ...
        'Color',            cmap(j,:), ...
        'MarkerFaceColor',  'none', ...
        'MarkerEdgeColor',  cmap(j,:), ...
        'MarkerSize',       9, ...
        'LineWidth',        1.8, ...
        'HandleVisibility', 'off');
end
lgd2 = legend(ax2, 'show', 'FontSize', 10, 'Box', 'off');
lgd2.Location = 'none';
lgd2.Position = [0.87 0.30 0.12 0.40];

%% ── FIGURE 3: a+g vs t ───────────────────────────────────────────────────
fig3 = figure('Name', 'a+g vs t', 'ToolBar', 'none', 'MenuBar', 'none');
fig3.Position = [100 100 980 520];
ax3  = axes(fig3, 'Position', [0.10 0.13 0.74 0.82]);
hold(ax3, 'on');

for j = 1:nH
    [~, marker] = get_height_style(heights(j).h_cm);
    for i = 1:heights(j).nTrials
        k        = heights(j).trials(i).kinematics;
        idx_plot = k.impact_index:5:k.stopFrame;
        plot(ax3, k.t_s(idx_plot), k.a_plus_g(idx_plot), marker, ...
            'Color',            [0 0 0 0.35], ...
            'MarkerSize',       4, ...
            'LineStyle',        'none', ...
            'HandleVisibility', 'off');
    end
end

for j = 1:nH
    hg  = heights(j);
    [tPhys, yMean, yStd, nValid] = interp_to_norm_time(hg, 'a_plus_g');
    lw  = 2.5; if hg.nTrials < 3, lw = 1.2; end
    plot_shaded_band(ax3, tPhys, yMean, yStd, nValid, cmap(j,:), ...
        sprintf('v_0 = %.0f \\pm %.0f cm s^{-1}', hg.v0_mean, hg.v0_std), 1, lw);
end

set(ax3, 'FontSize',   13,  'Box',        'on', ...
         'LineWidth',  1.2, 'XColor',     [0 0 0], ...
         'YColor',     [0 0 0], ...
         'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'TickDir',    'in', ...
         'XLim', [0 ax3.XLim(2)], 'YLim', [0 ax3.YLim(2)]);
grid(ax3, 'off');
xlabel(ax3, '$t$  (s)',              'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);
ylabel(ax3, '$a+g$  (cm s$^{-2}$)', 'FontSize', 16, 'Interpreter', 'latex', 'Color', [0 0 0]);
lgd3 = legend(ax3, 'show', 'FontSize', 10, 'Box', 'off');
lgd3.Location = 'none';
lgd3.Position = [0.84 0.28 0.14 0.50];

%% ── CLEAN LINE PLOTS ─────────────────────────────────────────────────────
fig5   = plot_mean_lines(heights, 'z_smooth',  '$z$  (cm)',              cmap);
fig6   = plot_mean_lines(heights, 'v_smooth',  '$v$  (cm s$^{-1}$)',    cmap, ...
             'yline', true, 'v0marks', true);
fig_ag = plot_mean_lines(heights, 'a_plus_g',  '$a+g$  (cm s$^{-2}$)', cmap, ...
             'xlim0', true, 'ylim0', true);

%% ── NORMALIZED COLLAPSE ──────────────────────────────────────────────────
fig7 = plot_normalized_collapse(heights, cmap);

%% ── a+g vs v0 ──────────────────────────────────────────────────
fig_dv0 = plot_depth_vs_v0(heights, cmap);

%% ── depth vs v0 ──────────────────────────────────────────────────
z_targets   = [0.5, 1.0, 1.5, 2.0];          % adjust to your depth range
cmap_depths = parula(numel(z_targets));
fig_3a      = plot_ag_vs_v2(heights, cmap_depths, z_targets);
%% ── SAVE ─────────────────────────────────────────────────────────────────
outDir = uigetdir(pwd, 'Select folder to save figures');
if outDir == 0
    warning('No folder selected — figures not saved.');
else
    fig1 = findobj('Type', 'figure', 'Name', 'z vs t');
    fig2 = findobj('Type', 'figure', 'Name', 'v vs t');

    figList  = {fig1, fig2, fig3, fig7, fig5, fig6, fig_ag};
    nameList = {'z_vs_t', 'v_vs_t', 'aplusg_vs_t', 'normalized_collapse', ...
                'z_vs_t_clean', 'v_vs_t_clean', 'aplusg_vs_t_clean'};

    for k = 1:numel(figList)
        base = fullfile(outDir, nameList{k});
        exportgraphics(figList{k}, [base '.pdf'], 'ContentType', 'vector');
        exportgraphics(figList{k}, [base '.png'], 'Resolution', 300);
    end

    save_analysis(heights, hVals, cmap, outDir);
end