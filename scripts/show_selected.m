function tags = show_selected(clearIt)
% SHOW_SELECTED  Print the curves clicked in dryrun_default_click, formatted
% ready to paste into a manual exclusion list.
%
%   show_selected        print the list
%   show_selected(true)  print it, then clear the selection

    if nargin < 1, clearIt = false; end
    tags = getappdata(0,'selectedExcludeTags');
    if isempty(tags)
        fprintf('No curves selected yet. Click curves in the dry-run figure.\n');
        tags = {}; return;
    end
    tags = unique(tags,'stable');

    fprintf('\n%d trial(s) selected for review:\n\n', numel(tags));
    fprintf('manualExclude = [ ...\n');
    for i = 1:numel(tags)
        sep = ', ...'; if i == numel(tags), sep = '];'; end
        fprintf('  "%s"%s\n', tags{i}, sep);
    end
    fprintf('\nplain list:\n');
    fprintf('%s\n', tags{:});

    if clearIt
        setappdata(0,'selectedExcludeTags',{});
        assignin('base','selectedExcludeTags',{});
        fprintf('\nselection cleared.\n');
    end
end
