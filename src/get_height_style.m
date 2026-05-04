function [col, marker, linestyle] = get_height_style(h_cm)

    styles = {
        15.24,  [0.20 0.10 0.60],  'o',  '-';    % circle
        20.32,  [0.15 0.35 0.80],  's',  '-';    % square
        25.40,  [0.10 0.60 0.70],  '^',  '-';    % triangle up
        30.48,  [0.15 0.70 0.30],  'd',  '-';    % diamond
        35.56,  [0.80 0.70 0.10],  'v',  '-.';   % triangle down
        40.64,  [0.85 0.40 0.10],  'p',  '-.';   % pentagon
        45.72,  [0.75 0.10 0.10],  'h',  '-.';   % hexagon
    };

    for k = 1:size(styles, 1)
        if abs(h_cm - styles{k,1}) < 1.5
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