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

%% Physical measurements: [h_cm, d_final_mm, d_std_mm]
physMeasured = [
%   h_cm    mean_mm   std_mm
    15.24,   15.5,     0.76;   
    20.32,   18.7,     0.33;
    25.40,   22.0,     0.70;
    30.48,   23.0,     0.55;
    35.56,   24.6,     0.69;
    40.64,   26.0,     0.82;
    45.72,   27.3,     1.24;
];
%% ── 2) FIGURE 1: z vs t ──────────────────────────────────────────────────
nH   = numel(heights);
cmap = cool(nH);   % sequential colormap: blue → magenta with height

[ax1, peakZ] = plot_time_series(heights, 'z_smooth', 'd_final_cm', 10, ...
    'z vs t', 'depth  (cm)', 'southeast', cmap);
set(ax1, 'YLim', [0, ax1.YLim(2)]);   % clip y at 0
add_peak_inset(peakZ, hVals, physMeasured, ...
    [0.55 0.15 0.32 0.32], 'h  (cm)', 'd_{final}  (mm)');

%% ── 3) FIGURE 2: v vs t ──────────────────────────────────────────────────
[ax2, peakV] = plot_time_series(heights, 'v_smooth', 'v0_cm_s', 1, ...
    'v vs t', 'v  (cm s^{-1})', 'northeast', cmap);
yline(ax2, 0, '--k', 'HandleVisibility', 'off');

% v0 markers at t=0 — hollow, per-height color from cmap
for j = 1:numel(heights)
    hg     = heights(j);
    col    = cmap(j,:);
    [~, marker] = get_height_style(hg.h_cm);
    errorbar(ax2, 0, hg.v0_mean, hg.v0_std, marker, ...
        'Color',           col, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', col, ...
        'MarkerSize',      9, ...
        'LineWidth',       1.8, ...
        'HandleVisibility','off');
end

% Inset: measured v0 (hollow colored) + theoretical sqrt(2gh) (hollow black)
ax_inset2 = add_peak_inset(peakV, hVals, [], ...
    [0.55 0.55 0.32 0.32], 'h  (cm)', 'v_0  (cm s^{-1})');
g_cm_s2 = 980;
for j = 1:numel(heights)
    hg  = heights(j);
    plot(ax_inset2, hg.h_cm, sqrt(2 * g_cm_s2 * hg.h_cm), 'd', ...
        'Color',           [0 0 0], ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerSize',      9, ...
        'LineWidth',       1.5, ...
        'HandleVisibility','off');
end
plot(ax_inset2, nan, nan, 'ks', 'MarkerFaceColor', 'none', ...
    'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'measured');
plot(ax_inset2, nan, nan, 'kd', 'MarkerFaceColor', 'none', ...
    'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', '\surd(2gh)');
legend(ax_inset2, 'show', 'FontSize', 8, 'Box', 'off', 'Location', 'northwest');

%% ── 4) FIGURE 3: a+g vs t ────────────────────────────────────────────────
fig3 = figure('Name', 'a+g vs t');
fig3.Position = [100 100 820 520];
ax3  = axes(fig3); hold(ax3, 'on');

% Raw dots
for j = 1:numel(heights)
    for i = 1:heights(j).nTrials
        k = heights(j).trials(i).kinematics;
        plot(ax3, k.t_s(k.impact_index+1:k.stopFrame), ...
             k.a_plus_g(k.impact_index+1:k.stopFrame), '.', ...
             'Color', [0 0 0 0.25], 'MarkerSize', 3, ...
             'HandleVisibility', 'off');
    end
end

% Shaded bands
peakAG = struct();
for j = 1:numel(heights)
    hg  = heights(j);
    col = cmap(j,:);
    [tPhys, yMean, yStd, nValid] = interp_to_norm_time(hg, 'a_plus_g');
    lw = 2.5; if hg.nTrials < 3, lw = 1.2; end
    plot_shaded_band(ax3, tPhys, yMean, yStd, nValid, col, ...
        sprintf('v_0 = %.0f \\pm %.0f cm s^{-1}', hg.v0_mean, hg.v0_std), 1, lw);
    peakAG.(sprintf('h%d', round(hg.h_cm))) = struct( ...
        'h',    hg.h_cm, ...
        'v0',   hg.v0_mean, ...
        'v0std',hg.v0_std, ...
        'vals', hg.d_final_cm * 10);   % placeholder — swap for ag_max if needed
end

xline(ax3, 0, '--k', 'HandleVisibility', 'off');
yline(ax3, 0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
legend(ax3, 'show', 'FontSize', 11, 'Box', 'off', 'Location', 'northeast');
set(ax3, 'FontSize', 13, 'Box', 'off', ...
    'XLim', [0, ax3.XLim(2)], ...   % clip x at 0
    'YLim', [0, ax3.YLim(2)]);      % clip y at 0
grid(ax3, 'on');
xlabel(ax3, 't  (s)',           'FontSize', 14);
ylabel(ax3, 'a+g  (cm s^{-2})', 'FontSize', 14);

% Inset: peak a+g vs v0
ax_inset3 = add_peak_inset_v0(peakAG, cmap, ...
    [0.55 0.55 0.32 0.32], 'v_0  (cm s^{-1})', '(a+g)_{max}  (cm s^{-2})');

%% ── 5) SAVE ──────────────────────────────────────────────────────────────
outDir = uigetdir(pwd, 'Select folder to save figures');
if outDir == 0
    warning('No folder selected — figures not saved.');
else
    save_all_figures({fig1, fig2, fig3}, ...
        {'z_vs_t', 'v_vs_t', 'aplusg_vs_t'}, outDir);
end
