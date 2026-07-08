function logPath = write_batch_log(rows, summ, logDir)
% WRITE_BATCH_LOG  Write a .txt summary of a batch run.
%
%   rows : struct array, one element per video, with fields
%          name, material, batch, dropHeight_mm, trialNum, container,
%          firstValidFrame, fps, status ('OK'|'FAILED'), reason,
%          framesDir, detDir, resultsDir
%   summ : struct with total, nSuccess, nFailed, inputDir, timestamp
%   logDir : folder to write the log into

    if ~exist(logDir, 'dir'), mkdir(logDir); end
    stamp   = datestr(now, 'yyyymmdd_HHMMSS');
    logPath = fullfile(logDir, sprintf('batch_log_%s.txt', stamp));
    fid = fopen(logPath, 'w');

    fprintf(fid, '============================================================\n');
    fprintf(fid, ' BATCH PROCESSING LOG\n');
    fprintf(fid, ' Run time      : %s\n', summ.timestamp);
    fprintf(fid, ' Input folder  : %s\n', summ.inputDir);
    fprintf(fid, '------------------------------------------------------------\n');
    fprintf(fid, ' Total videos found : %d\n', summ.total);
    fprintf(fid, ' Successful         : %d\n', summ.nSuccess);
    fprintf(fid, ' Failed / skipped   : %d\n', summ.nFailed);
    fprintf(fid, '============================================================\n\n');

    for i = 1:numel(rows)
        r = rows(i);
        fprintf(fid, '[%3d] %-26s  %s\n', i, r.name, r.status);
        fprintf(fid, '      material=%s  batch=%s  drop=%gmm  trial=T%02d  container=%s\n', ...
            r.material, r.batch, r.dropHeight_mm, r.trialNum, r.container);
        fprintf(fid, '      firstValidFrame=%s  fps=%s\n', ...
            num2str(r.firstValidFrame), fmtnum(r.fps, '%.4f'));
        if strcmpi(r.status, 'OK')
            fprintf(fid, '      frames  : %s\n', r.framesDir);
            fprintf(fid, '      detect  : %s\n', r.detDir);
            fprintf(fid, '      results : %s\n', r.resultsDir);
        end
        if ~isempty(r.reason)
            fprintf(fid, '      note    : %s\n', r.reason);
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('\nBatch log written: %s\n', logPath);
end

function s = fmtnum(x, fmt)
    if isnan(x), s = 'NaN'; else, s = sprintf(fmt, x); end
end
