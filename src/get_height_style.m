function [col, marker, linestyle] = get_height_style(h_cm)
% GET_HEIGHT_STYLE  Return color, marker, and linestyle for a given drop height.
%
%   Supports any number of heights by mapping h_cm to a fixed set of
%   markers and linestyles via sorted rank. Color is set to a parula-based
%   value but is typically overridden by the cmap colorbar in overlay plots.
%
%   Tolerance: ±1.5 cm for height matching.

    % ── All recognised heights (cm) ───────────────────────────────────────
    knownH = [7.62, 10.16, 12.70, 15.24, 17.78, 20.32, 22.86, ...
              25.40, 27.94, 30.48, 33.02, 35.56, 38.10, 40.64, 43.18, 45.72];

    % ── Marker cycle (enough for 16 heights without repeating) ────────────
    markers    = {'o','s','d','^','v','p','h','>','<','*', ...
                  'o','s','d','^','v','p'};
    linestyles = {'-','-','-','-','-','-','-','-', ...
                  '-.','-.','-.','-.','-.','-.','-.', '-.'};

    % ── Match ─────────────────────────────────────────────────────────────
    [minDist, idx] = min(abs(h_cm - knownH));

    if minDist > 1.5
        % Truly unrecognised — grey fallback
        col       = [0.40 0.40 0.40];
        marker    = 'x';
        linestyle = ':';
        warning('get_height_style: unrecognized h=%.2f cm — using grey fallback.', h_cm);
        return;
    end

    % ── Color: parula-based, evenly spaced across known heights ───────────
    nH  = numel(knownH);
    cm  = parula(nH);
    col = cm(idx, :);

    marker    = markers{idx};
    linestyle = linestyles{idx};
end