function L = batch_log_init(logDir, stamp, prefixes, headerText, csvHeader)
% BATCH_LOG_INIT  Shared crash-safe batch logger — open txt+csv, write headers.
%   prefixes = {txtPrefix, csvPrefix, retryPrefix}; filenames get '<stamp>.<ext>'.
%   headerText/csvHeader are pre-formatted (real newlines) by the caller, so the
%   column layout stays stage-specific while the file mechanics are shared.
    if ~exist(logDir, 'dir'), mkdir(logDir); end
    L.stamp = stamp;
    L.txt   = fullfile(logDir, sprintf('%s%s.txt', prefixes{1}, stamp));
    L.csv   = fullfile(logDir, sprintf('%s%s.csv', prefixes{2}, stamp));
    L.retry = fullfile(logDir, sprintf('%s%s.txt', prefixes{3}, stamp));
    fid = fopen(L.txt, 'w'); fprintf(fid, '%s', headerText); fclose(fid);
    fid = fopen(L.csv, 'w'); fprintf(fid, '%s', csvHeader);  fclose(fid);
end
