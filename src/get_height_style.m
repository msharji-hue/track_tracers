function [col, marker, linestyle] = get_height_style(h_cm)
    if abs(h_cm - 15.24) < 1
        col       = [0.00 0.27 0.60];   % deep blue
        marker    = 'o';
        linestyle = '-';
    elseif abs(h_cm - 30.48) < 1
        col       = [0.85 0.45 0.00];   % orange
        marker    = 's';
        linestyle = '-';
    elseif abs(h_cm - 45.72) < 1
        col       = [0.70 0.09 0.12];   % crimson
        marker    = 'd';
        linestyle = '-.';
    else
        col       = [0.40 0.40 0.40];
        marker    = 'x';
        linestyle = ':';
    end
end