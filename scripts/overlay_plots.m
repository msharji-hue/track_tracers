%% overlay_plots.m
clear; close all; clc;

codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

%% 1) LOAD TRIALS
trials  = load_selected_trials();
hVals   = unique([trials.h_cm]);
nTrials = numel(trials);

physMeasured = [15.24, 13.0, 1.30;
                30.48, 21.2, 0.96;
                45.72, 25.0, 1.67];

d1_cm = 8.7;

%% 1.5) SELECT ONE REPRESENTATIVE TRIAL PER DROP HEIGHT
repTrials = trials(1:0);   % empty struct array with same fields as trials
for j = 1:numel(hVals)
    h   = hVals(j);
    idx = find([trials.h_cm] == h);
    fprintf('\nTrials for h = %.2f cm:\n', h);
    for k = 1:numel(idx)
        fprintf('  [%d] T%02d  v0=%.1f cm/s  d=%.3f cm\n', k, ...
            trials(idx(k)).trialInfo.trialNum, ...
            trials(idx(k)).scalars.v0_cm_s, ...
            trials(idx(k)).scalars.d_final_cm);
    end
    choice = input(sprintf('Select trial for overlay (1-%d): ', numel(idx)));
    repTrials(end+1) = trials(idx(choice));
end
nOverlay = numel(repTrials);
%% ── FIGURE 1: z vs t ────────────────────────────────────────────────────
fig1 = figure('Name', 'z vs t');
ax1  = axes(fig1); hold(ax1, 'on');
peakZ = plot_overlay_lines(ax1, repTrials, trials, hVals, 'z_smooth', 'd_final_cm', 10);
xline(ax1, 0, '--k', 'impact', 'HandleVisibility', 'off');
legend(ax1, 'show', 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');
set(ax1, 'FontSize', 13, 'Box', 'off'); grid(ax1, 'on');
xlabel(ax1, 't (s)', 'FontSize', 14); ylabel(ax1, 'depth (cm)', 'FontSize', 14);
add_peak_inset(peakZ, hVals, physMeasured, [0.55 0.15 0.32 0.32], 'h (cm)', 'd_{final} (mm)');

%% ── FIGURE 2: v vs t ────────────────────────────────────────────────────
fig2 = figure('Name', 'v vs t');
ax2  = axes(fig2); hold(ax2, 'on');
peakV = plot_overlay_lines(ax2, repTrials, trials, hVals, 'v_smooth', 'v0_cm_s', 1);
xline(ax2, 0, '--k', 'HandleVisibility', 'off');
yline(ax2, 0, '--k', 'HandleVisibility', 'off');
legend(ax2, 'show', 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');
set(ax2, 'FontSize', 13, 'Box', 'off'); grid(ax2, 'on');
xlabel(ax2, 't (s)', 'FontSize', 14); ylabel(ax2, 'v (cm/s)', 'FontSize', 14);
add_peak_inset(peakV, hVals, [], [0.55 0.55 0.32 0.32], 'h (cm)', 'v_0 (cm/s)');

%% ── FIGURE 3: a+g vs t ──────────────────────────────────────────────────
fig3 = figure('Name', 'a+g vs t');
ax3  = axes(fig3); hold(ax3, 'on');
plot_overlay_lines(ax3, repTrials, trials, hVals, 'a_plus_g', 'v0_cm_s', 1);
xline(ax3, 0, '--k', 'HandleVisibility', 'off');
yline(ax3, 0, '--k', 'HandleVisibility', 'off');
legend(ax3, 'show', 'FontSize', 12, 'Box', 'off');
set(ax3, 'FontSize', 13, 'Box', 'off'); grid(ax3, 'on');
xlabel(ax3, 't (s)', 'FontSize', 14); ylabel(ax3, 'a+g (cm/s²)', 'FontSize', 14);
%% ── FIGURE 4: v vs depth ─────────────────────────────────────────────────
fig4 = figure('Name', 'v vs depth');
ax4  = axes(fig4, 'Position', [0.10 0.11 0.65 0.82]);
hold(ax4, 'on');

for i = 1:nOverlay
    k      = repTrials(i).kinematics;
    h      = repTrials(i).h_cm;
    [~, marker] = get_height_style(h);
    t_norm = (k.t_s - k.t_s(k.impact_index)) ./ ...
             (k.t_s(k.stopFrame) - k.t_s(k.impact_index));
    t_norm = max(0, min(1, t_norm));
    scatter(ax4, k.z_smooth, k.v_smooth, 22, t_norm, marker, 'filled', ...
        'HandleVisibility', 'off');
end

cb = colorbar(ax4); colormap(ax4, jet);
cb.Label.String = 'normalized time';
cb.Label.FontSize = 12;

for h = hVals
    [~, marker] = get_height_style(h);
    scatter(ax4, nan, nan, 50, [0.4 0.4 0.4], marker, 'filled', ...
        'DisplayName', sprintf('h = %.2f cm', h));
end
legend(ax4, 'show', 'FontSize', 12, 'Box', 'off', 'Location', 'northeast');
set(ax4, 'FontSize', 13, 'Box', 'off'); grid(ax4, 'on');
xlabel(ax4, 'depth (cm)', 'FontSize', 14);
ylabel(ax4, 'v (cm/s)',   'FontSize', 14);

%% ── FIGURE 5: normalized collapse v/v0 vs z/z_max ───────────────────────
% Compute alpha fit + pre-store per-trial normalized data
z_all    = [];
v_all    = [];
trial_z  = cell(nOverlay, 1);
trial_v  = cell(nOverlay, 1);

for i = 1:nOverlay
    k     = repTrials(i).kinematics;
    v0    = k.v_smooth(k.impact_index);
    z_max = max(k.z_smooth, [], 'omitnan');
    if v0 == 0 || z_max == 0; continue; end
    z_n   = k.z_smooth ./ z_max;
    v_n   = k.v_smooth ./ v0;
    valid = isfinite(z_n) & isfinite(v_n) & z_n >= 0 & z_n <= 1;
    trial_z{i} = z_n(valid);
    trial_v{i} = v_n(valid);
    z_all = [z_all; z_n(valid)];
    v_all = [v_all; v_n(valid)];
end

mask      = z_all < 0.99 & v_all > 0;
alpha_fit = log(1 - z_all(mask)) \ log(v_all(mask));

% Trial-level bootstrap 90% CI (honest: only 3 trials)
n_boot     = 500;
alpha_boot = zeros(n_boot, 1);
rng(42);
for kb = 1:n_boot
    idx_t  = randi(nOverlay, nOverlay, 1);       % resample trials w/ replacement
    z_s    = vertcat(trial_z{idx_t});
    v_s    = vertcat(trial_v{idx_t});
    m      = z_s < 0.99 & v_s > 0;
    alpha_boot(kb) = log(1 - z_s(m)) \ log(v_s(m));
end
alpha_lo = prctile(alpha_boot, 5);
alpha_hi = prctile(alpha_boot, 95);

fig5 = figure('Name', 'normalized collapse');
fig5.Position = [100 100 800 550];
ax5  = axes(fig5); hold(ax5, 'on');

% Bootstrap confidence band (draw first so data plots on top)
zRef = linspace(0, 1, 300);
v_lo = max(0, 1 - zRef).^alpha_hi;
v_hi = max(0, 1 - zRef).^alpha_lo;
fill(ax5, [zRef, fliplr(zRef)], [v_hi, fliplr(v_lo)], ...
    [0.3 0.3 0.3], 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

% Data
for i = 1:nOverlay
    k     = repTrials(i).kinematics;
    h     = repTrials(i).h_cm;
    col   = get_height_style(h);
    v0    = k.v_smooth(k.impact_index);
    z_max = max(k.z_smooth, [], 'omitnan');
    if v0 == 0 || z_max == 0; continue; end
    plot(ax5, k.z_smooth ./ z_max, k.v_smooth ./ v0, '-', ...
        'Color', col, 'LineWidth', 2.5, 'HandleVisibility', 'off');
end

% Reference curves
plot(ax5, zRef, sqrt(max(0, 1 - zRef)), 'k:', ...
    'LineWidth', 1.5, 'DisplayName', '\surd(1-z/z_{max})  [const. drag]');
plot(ax5, zRef, max(0, 1 - zRef).^alpha_fit, '--', ...
    'Color', [0.15 0.15 0.15], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('(1-z/z_{max})^{%.2f}  [fit]', alpha_fit));

% Height legend entries
for h = hVals
    col = get_height_style(h);
    plot(ax5, nan, nan, '-', 'Color', col, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('h = %.2f cm', h));
end

legend(ax5, 'show', 'FontSize', 12, 'Box', 'off', 'Location', 'northeast');
set(ax5, 'FontSize', 13, 'Box', 'off', 'XLim', [0 1], 'YLim', [0 1.05]);
xlabel(ax5, 'z / z_{max}', 'FontSize', 14);
ylabel(ax5, 'v / v_0',     'FontSize', 14);
grid(ax5, 'on');

text(ax5, 0.97, 0.50, ...
    sprintf('\\alpha = %.2f  [%.2f, %.2f]', alpha_fit, alpha_lo, alpha_hi), ...
    'FontSize', 12, 'Units', 'normalized', 'FontAngle', 'italic', ...
    'Color', [0.2 0.2 0.2], 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');

%% ── FIGURE 6: drag phase portrait — colored by normalized depth ──────────
fig6 = figure('Name', 'drag phase portrait');
fig6.Position = [100 100 850 600];
ax6  = axes(fig6, 'Position', [0.10 0.11 0.70 0.82]); hold(ax6, 'on');

for i = 1:nOverlay
    k     = repTrials(i).kinematics;
    idx   = k.impact_index : k.stopFrame;
    v     = k.v_smooth(idx);
    ag    = k.a_plus_g(idx);
    z     = k.z_smooth(idx);
    z_max = max(z, [], 'omitnan');
    valid = isfinite(v) & isfinite(ag) & isfinite(z) & v >= 0 & z > 0 & ag > 0;
    z_n   = z(valid) ./ z_max;
    scatter(ax6, v(valid), ag(valid), 35, z_n, 'filled', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 0.3, 'HandleVisibility', 'off');
end

colormap(ax6, turbo);
cb = colorbar(ax6);
cb.Label.String   = 'z / z_{max}';
cb.Label.FontSize = 13;
cb.FontSize       = 11;
clim(ax6, [0 1]);

yline(ax6, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
    'HandleVisibility', 'off');
set(ax6, 'FontSize', 13, 'Box', 'off', ...
    'YLim', [0, max(ax6.YLim)], ...
    'XLim', [0 300]);
grid(ax6, 'on');
xlabel(ax6, 'v  (cm s^{-1})',     'FontSize', 14);
ylabel(ax6, 'a + g  (cm s^{-2})', 'FontSize', 14);
text(ax6, 0.97, 0.05, ...
    {'At equal v, force increases with depth', ...
     '\Rightarrow  F \neq f(v) alone'}, ...
    'FontSize', 11, 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'Color', [0.25 0.25 0.25], 'FontAngle', 'italic');
%% ── SAVE ─────────────────────────────────────────────────────────────────
outDir = uigetdir(pwd, 'Select folder to save overlay figures');
if outDir == 0
    warning('No folder selected — figures not saved.');
else
    figs  = {fig1, fig2, fig3, fig4, fig5, fig6};
    names = {'z_vs_t', 'v_vs_t', 'aplusg_vs_t', 'v_vs_depth', ...
             'normalized_collapse', 'drag_phase_portrait'};

    % ── sub-folders ───────────────────────────────────────────────────────
    pngDir = fullfile(outDir, 'png');
    pdfDir = fullfile(outDir, 'pdf');
    figDir = fullfile(outDir, 'fig');
    mkdir(pngDir); mkdir(pdfDir); mkdir(figDir);

    for f = 1:numel(figs)
        base = names{f};
        exportgraphics(figs{f}, fullfile(pngDir, [base '.png']), ...
            'Resolution', 300, 'BackgroundColor', 'white');
        exportgraphics(figs{f}, fullfile(pdfDir, [base '.pdf']), ...
            'ContentType', 'vector',  'BackgroundColor', 'white');
        savefig(figs{f}, fullfile(figDir, [base '.fig']));
        fprintf('Saved: %s\n', base);
    end
    fprintf('\nAll figures saved to:\n  %s\n', outDir);
end
