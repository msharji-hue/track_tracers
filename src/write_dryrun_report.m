function write_dryrun_report(outRoot, inputDesc, rows)
% WRITE_DRYRUN_REPORT  Shared dry-run reporter (console + timestamped .txt).
%   rows : struct array with fields
%            .head       char, the one-line summary for this item
%            .ok         logical
%            .pathLines  cell{N x 2} of {label,path}, printed only when ok
%   Reproduces the process_trial dry-run format exactly; used by both stages.
    stamp  = datestr(now, 'yyyymmdd_HHMMSS');
    repDir = fullfile(outRoot, '03_RESULTS', '_batch_logs');
    if ~exist(repDir, 'dir'), mkdir(repDir); end
    repPath = fullfile(repDir, sprintf('dryrun_report_%s.txt', stamp));
    fid = fopen(repPath, 'w');

    n = numel(rows);
    fprintf('DRY RUN — %d videos found under %s\n\n', n, inputDesc);
    fprintf(fid, 'DRY RUN  %s\nInput: %s\nVideos found: %d\n\n', stamp, inputDesc, n);

    nOK = 0; nBad = 0;
    for i = 1:n
        r = rows(i);
        fprintf('%s\n', r.head);
        if r.ok
            nOK = nOK + 1;
            fprintf(fid, '%s\n', r.head);
            for k = 1:size(r.pathLines, 1)
                fprintf(fid, '      %-7s: %s\n', r.pathLines{k,1}, r.pathLines{k,2});
            end
            fprintf(fid, '\n');
        else
            nBad = nBad + 1;
            fprintf(fid, '%s\n\n', r.head);
        end
    end
    fprintf('\nDry run: %d processable, %d unparseable.\nReport: %s\n', nOK, nBad, repPath);
    fprintf(fid, 'SUMMARY: %d processable, %d unparseable\n', nOK, nBad);
    fclose(fid);
end
