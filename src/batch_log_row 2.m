function batch_log_row(L, txtBlock, csvRow, retryLine)
% BATCH_LOG_ROW  Append one trial's pre-formatted txt block + csv row, closing
%   each file every call (crash-safe). Optional retryLine appended to L.retry.
    fid = fopen(L.txt, 'a'); fprintf(fid, '%s', txtBlock); fclose(fid);
    fid = fopen(L.csv, 'a'); fprintf(fid, '%s', csvRow);   fclose(fid);
    if nargin >= 4 && ~isempty(retryLine)
        fid = fopen(L.retry, 'a'); fprintf(fid, '%s\n', retryLine); fclose(fid);
    end
end
