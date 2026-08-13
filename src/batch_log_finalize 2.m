function batch_log_finalize(L, txtBlock)
% BATCH_LOG_FINALIZE  Append the pre-formatted summary block to the txt log.
    fid = fopen(L.txt, 'a'); fprintf(fid, '%s', txtBlock); fclose(fid);
end
