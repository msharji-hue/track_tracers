function plot_drag_portrait(ax, repTrials, nOverlay, hVals)
% PLOT_DRAG_PORTRAIT  Scatter of a+g vs v, colored by normalized depth, shape by height.

    vMax  = max(cellfun(@(t) max(t.kinematics.v_smooth(t.kinematics.impact_index:t.kinematics.stopFrame),  [], 'omitnan'), num2cell(repTrials)));
    agMax = max(cellfun(@(t) max(t.kinematics.a_plus_g(t.kinematics.impact_index:t.kinematics.stopFrame), [], 'omitnan'), num2cell(repTrials)));

    for i = 1:nOverlay
        k     = repTrials(i).kinematics;
        [~, marker] = get_height_style(repTrials(i).h_cm);
        idx   = k.impact_index:k.stopFrame;
        v     = k.v_smooth(idx);  ag = k.a_plus_g(idx);  z = k.z_smooth(idx);
        z_max = max(z, [], 'omitnan');
        valid = isfinite(v) & isfinite(ag) & isfinite(z) & v >= 0 & z > 0 & ag > 0;
        scatter(ax, v(valid), ag(valid), 35, z(valid)./z_max, marker, 'filled', ...
            'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
            'LineWidth', 0.3, 'HandleVisibility', 'off');
    end

    colormap(ax, turbo);
    cb = colorbar(ax); cb.Label.String = 'z / z_{max}'; cb.Label.FontSize = 13; cb.FontSize = 11;
    clim(ax, [0 1]);

    for h = hVals
        [~, marker] = get_height_style(h);
        scatter(ax, nan, nan, 50, [0.5 0.5 0.5], marker, 'filled', 'DisplayName', sprintf('h = %.2f cm', h));
    end
    legend(ax, 'show', 'FontSize', 11, 'Box', 'off', 'Location', 'northwest');
    set(ax, 'FontSize', 13, 'Box', 'off', 'XLim', [0 vMax*1.05], 'YLim', [0 agMax*1.05]);
    grid(ax, 'on');
    xlabel(ax, 'v  (cm s^{-1})', 'FontSize', 14); ylabel(ax, 'a + g  (cm s^{-2})', 'FontSize', 14);
    text(ax, 0.97, 0.05, {'At equal v, force increases with depth', '\Rightarrow  F \neq f(v) alone'}, ...
        'FontSize', 11, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontAngle', 'italic');
end
