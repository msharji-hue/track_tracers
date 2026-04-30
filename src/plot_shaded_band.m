function plot_shaded_band(ax, tPhys, yMean, yStd, nValid, col, dispName, minN, lw)
% PLOT_SHADED_BAND  Mean line + ±1 std shaded band.
%
%   Inputs:
%       ax       - axes handle
%       tPhys    - physical time vector (s)
%       yMean    - mean values
%       yStd     - std values
%       nValid   - number of valid trials per timepoint
%       col      - [R G B] color
%       dispName - legend label string
%       minN     - minimum trials to plot (default 3)
%       lw       - line width (default 2.5, use 1.2 for sparse data)

    if nargin < 8, minN = 3;   end
    if nargin < 9, lw   = 2.5; end

    ok = nValid >= minN & isfinite(yMean) & isfinite(yStd);
    if ~any(ok), warning('plot_shaded_band: no valid data to plot.'); return; end

    t_ok = tPhys(ok);
    m_ok = yMean(ok);
    s_ok = yStd(ok);

    % Shaded ±1 std band
    fill(ax, [t_ok, fliplr(t_ok)], ...
         [m_ok + s_ok, fliplr(m_ok - s_ok)], ...
         col, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');

    % Mean line
    plot(ax, t_ok, m_ok, '-', ...
        'Color',       col, ...
        'LineWidth',   lw, ...
        'DisplayName', dispName);
end