function apply_fig_style(ax)
%APPLY_FIG_STYLE  The shared axes convention for this repo's publication figures.
%
%   apply_fig_style(ax)
%
%   Applies exactly what scripts/depth_scaling.m has always used, so the
%   kinematics figures and the scaling figures stay visually matched:
%
%       hold on, grid on, box on
%
%   and NOTHING ELSE. Tick direction, minor ticks, tick length, font size and
%   font name are left at MATLAB's defaults deliberately -- that is what
%   depth_scaling.m does, and overriding any of them here would silently
%   desynchronise the two families of figure the moment one script was edited
%   and the other was not.
%
%   Tick LABELS are left alone too: callers that want to blank them must do so
%   themselves. The kinematics figures deliberately keep full numeric labels on
%   both axes of every panel.
%
%   If the house style ever changes, change it here and both scripts follow.

    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');
end
