%% scripts/overlay_plots.m
% Overlay plots: z vs t, v vs t, a+g vs t
% Shaded band = mean ± 1 std across ~5 trials per height
% X-axis in physical time (s), legend shows v0 not h
clear; close all; clc;

codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

%% LOAD + GROUP
trials  = load_selected_trials();
heights = group_trials_by_height(trials);
hVals   = [heights.h_cm];

%% ── COLORMAP ─────────────────────────────────────────────────────────────
nH         = numel(heights);
cmap_clean = parula(nH);
cmap       = cmap_clean;   % use same colormap for all figures

%% ── 2) FIGURE 1: z vs t ──────────────────────────────────────────────────
[ax1, peakZ] = plot_time_series(heights, 'z_smooth', 'd_final_cm', 10, ...
    'z vs t', 'depth  (cm)', 'none', cmap);

set(ax1, 'Box', 'on', 'LineWidth', 1.2, 'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);

lgd = legend(ax1, 'show', 'FontSize', 10, 'Box', 'off');
lgd.Location = 'none';
lgd.Position = [0.87 0.30 0.12 0.40];

ax_inset1 = add_peak_inset(peakZ, [0.55 0.22 0.22 0.26], ...
    'v_0  (cm s^{-1})', 'd_{final}  (mm)');
set(ax_inset1, 'LineWidth', 1.2, 'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1], ...
    'YLim', [12, ax_inset1.YLim(2)]);

%% ── 3) FIGURE 2: v vs t ──────────────────────────────────────────────────
[ax2, peakV] = plot_time_series(heights, 'v_smooth', 'v0_cm_s', 1, ...
    'v vs t', 'v  (cm s^{-1})', 'none', cmap);
yline(ax2, 0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'HandleVisibility', 'off');

for j = 1:numel(heights)
    hg          = heights(j);
    col         = cmap(j,:);
    [~, marker] = get_height_style(hg.h_cm);
    errorbar(ax2, 0, hg.v0_mean, hg.v0_std, marker, ...
        'Color',           col, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', col, ...
        'MarkerSize',      9, ...
        'LineWidth',       1.8, ...
        'HandleVisibility','off');
end

ax_inset2 = axes('Position', [0.40 0.50 0.26 0.30]);
hold(ax_inset2, 'on');

for j = 1:numel(heights)
    hg          = heights(j);
    col         = cmap(j,:);
    [~, marker] = get_height_style(hg.h_cm);
    tstops = arrayfun(@(t) t.scalars.t_stop_s, hg.trials) * 1000;
    v0s    = arrayfun(@(t) t.scalars.v0_cm_s,  hg.trials);
    errorbar(ax_inset2, mean(v0s), mean(tstops), std(tstops), std(tstops), ...
        std(v0s), std(v0s), marker, ...
        'Color',           col, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', col, ...
        'MarkerSize',      8, ...
        'LineWidth',       1.8, ...
        'HandleVisibility','off');
end

set(ax_inset2, 'FontSize', 10, 'Box', 'on', 'LineWidth', 1.2, ...
    'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);
xlabel(ax_inset2, 'v_0  (cm s^{-1})', 'FontSize', 11);
ylabel(ax_inset2, 't_{stop}  (ms)',    'FontSize', 11);
grid(ax_inset2, 'off');

%% ── 4) FIGURE 3: a+g vs t ────────────────────────────────────────────────
fig3 = figure('Name', 'a+g vs t');
fig3.Position = [100 100 980 520];
ax3  = axes(fig3, 'Position', [0.08 0.11 0.74 0.82]); hold(ax3, 'on');

for j = 1:numel(heights)
    [~, marker] = get_height_style(heights(j).h_cm);
    for i = 1:heights(j).nTrials
        k        = heights(j).trials(i).kinematics;
        idx_plot = k.impact_index:5:k.stopFrame;
        plot(ax3, k.t_s(idx_plot), ...
             k.a_plus_g(idx_plot), marker, ...
             'Color',            [0 0 0 0.35], ...
             'MarkerSize',       4, ...
             'LineStyle',        'none', ...
             'HandleVisibility', 'off');
    end
end

peakAG = struct();
for j = 1:numel(heights)
    hg  = heights(j);
    col = cmap(j,:);
    [tPhys, yMean, yStd, nValid] = interp_to_norm_time(hg, 'a_plus_g');
    lw = 2.5; if hg.nTrials < 3, lw = 1.2; end
    plot_shaded_band(ax3, tPhys, yMean, yStd, nValid, col, ...
        sprintf('v_0 = %.0f \\pm %.0f cm s^{-1}', hg.v0_mean, hg.v0_std), 1, lw);
    ag_peaks = arrayfun(@(t) max(t.kinematics.a_plus_g( ...
        t.kinematics.impact_index:t.kinematics.stopFrame), [], 'omitnan'), hg.trials);
    peakAG.(sprintf('h%d', round(hg.h_cm))) = struct( ...
        'h', hg.h_cm, 'v0', hg.v0_mean, 'v0std', hg.v0_std, 'vals', ag_peaks);
end

lgd          = legend(ax3, 'show', 'FontSize', 10, 'Box', 'off');
lgd.Location = 'none';
lgd.Position = [0.84 0.28 0.14 0.50];

set(ax3, 'FontSize', 13, 'Box', 'on', 'LineWidth', 1.2, ...
    'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1], ...
    'XLim',   [0, ax3.XLim(2)], 'YLim', [0, ax3.YLim(2)]);
grid(ax3, 'off');
xlabel(ax3, 't  (s)',           'FontSize', 14);
ylabel(ax3, 'a+g  (cm s^{-2})', 'FontSize', 14);

ax_inset3 = add_peak_inset(peakAG, [0.50 0.50 0.30 0.30], ...
    'v_0  (cm s^{-1})', '(a+g)_{max}  (cm s^{-2})');
set(ax_inset3, 'LineWidth', 1.2, 'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);

%% ── 5) SAVE ──────────────────────────────────────────────────────────────
outDir = uigetdir(pwd, 'Select folder to save figures');
if outDir == 0
    warning('No folder selected — figures not saved.');
else
    % Get figure handles by name
    fig1 = findobj('Type', 'figure', 'Name', 'z vs t');
    fig2 = findobj('Type', 'figure', 'Name', 'v vs t');
    % fig3 already defined

    save_all_figures({fig1, fig2, fig3}, ...
        {'z_vs_t', 'v_vs_t', 'aplusg_vs_t'}, outDir);
    save_analysis(heights, hVals, cmap, peakZ, peakV, peakAG, outDir);
end