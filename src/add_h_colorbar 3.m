function cb = add_h_colorbar(ax, hVals, cmap)
% ADD_H_COLORBAR  Add a discrete colorbar keyed to drop height values.
%
%   Inputs:
%       ax    - axes handle to attach colorbar to
%       hVals - vector of h_cm values (one per height group)
%       cmap  - [nH x 3] colormap matching hVals order

    nH = numel(hVals);
    colormap(ax, cmap);
    clim(ax, [1 nH]);
    cb = colorbar(ax, 'Location', 'eastoutside');
    cb.Ticks      = 1:nH;
    cb.TickLabels = arrayfun(@(h) sprintf('%.0f', h), hVals, 'UniformOutput', false);
    cb.Label.String   = 'h  (cm)';
    cb.Label.FontSize = 12;
    cb.FontSize       = 10;
end