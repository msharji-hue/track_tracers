function [col, marker, linestyle] = get_height_style(h_cm)

    styles = {
        15.24,  [0.00 0.27 0.60],  '^',  '-';    % deep blue
        20.32,  [0.00 0.65 0.65],  'o',  '-';    % teal
        25.40,  [0.20 0.70 0.30],  's',  '-';    % green
        30.48,  [0.85 0.45 0.00],  'v',  '-';    % orange
        35.56,  [0.75 0.10 0.50],  'x',  '-.';   % purple
        40.64,  [0.85 0.10 0.80],  'd',  '-.';   % magenta
        45.72,  [0.70 0.09 0.12],  'p',  '-.';   % crimson
    };

    for k = 1:size(styles, 1)
        if abs(h_cm - styles{k,1}) < 1.5    % <-- widened tolerance to 1.5 cm
            col       = styles{k,2};
            marker    = styles{k,3};
            linestyle = styles{k,4};
            return;
        end
    end

    col       = [0.40 0.40 0.40];
    marker    = 'x';
    linestyle = ':';
    warning('get_height_style: unrecognized h=%.2f cm — using grey fallback.', h_cm);
end