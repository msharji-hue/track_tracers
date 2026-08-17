function save_tracks(meta, tracks, calib, outRoot)
% SAVE_TRACKS  Save tracked marker positions + metadata (tracking-only).
%
%   Writes into  <outRoot>/tracks/ :
%       <trialTag>_tracks.mat    meta + tracks + calib
%       <trialTag>_tracks.csv    per-frame x/y for every marker (wide format)
%       <trialTag>_scalars.csv   one-row trial metadata summary
%
%   No kinematics here: velocity, acceleration, (a+g), event indices, and rod
%   bending are computed DOWNSTREAM from these tracks.
%
%   Indexing: trackedX/trackedY columns use TRACKING-frame numbering, where
%   column 1 == firstValidFrame in the original video. The per-frame CSV
%   carries both 'tracking_frame' and 'original_frame'.

    trackDir = fullfile(outRoot, 'tracks');
    if ~exist(trackDir, 'dir'), mkdir(trackDir); end
    base = meta.trialTag;

    % ── Full bundle (.mat) ────────────────────────────────────────────────
    outMat = fullfile(trackDir, [base '_tracks.mat']);
    save(outMat, 'meta', 'tracks', 'calib', '-v7.3');
    fprintf('  Saved tracks MAT : %s\n', outMat);

    % ── Per-frame positions (.csv, wide) ─────────────────────────────────
    outCsvT = fullfile(trackDir, [base '_tracks.csv']);
    nM = size(tracks.trackedX, 1);
    nF = size(tracks.trackedX, 2);
    fid = fopen(outCsvT, 'w');

    hdr = 'tracking_frame,original_frame';
    for mm = 1:nM, hdr = [hdr sprintf(',x%d', mm)]; end %#ok<AGROW>
    for mm = 1:nM, hdr = [hdr sprintf(',y%d', mm)]; end %#ok<AGROW>
    fprintf(fid, '%s\n', hdr);

    % firstValidFrame indexes the EXPORTED window, so the absolute video frame
    % needs the window offset too. windowStart is 1 for a full-range export, in
    % which case this reduces to the old firstValidFrame + f - 1.
    wStart = 1;
    if isfield(meta,'windowStart') && isfinite(meta.windowStart)
        wStart = meta.windowStart;
    end
    for f = 1:nF
        origF = wStart + meta.firstValidFrame + f - 2;
        row   = sprintf('%d,%d', f, origF);
        for mm = 1:nM, row = [row sprintf(',%.4f', tracks.trackedX(mm,f))]; end %#ok<AGROW>
        for mm = 1:nM, row = [row sprintf(',%.4f', tracks.trackedY(mm,f))]; end %#ok<AGROW>
        fprintf(fid, '%s\n', row);
    end
    fclose(fid);
    fprintf('  Saved tracks CSV : %s\n', outCsvT);

    % ── One-row trial metadata (.csv) ────────────────────────────────────
    outCsv = fullfile(trackDir, [base '_scalars.csv']);
    fid = fopen(outCsv, 'w');
    fprintf(fid, ['material,batch,dropHeight_mm,trialNum,condition,', ...
                  'rho_particle_g_cm3,rho_bulk_g_cm3,phi,', ...
                  'fps,nFrames,firstValidFrame,nTracked,nMarkers,', ...
                  'windowStart,windowEnd,autoWindow\n']);
    fprintf(fid, '%s,%s,%g,%d,%s,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%g,%g,%d\n', ...
        meta.material, meta.batchName, meta.dropHeight_mm, meta.trialNum, ...
        meta.container, fld(meta,'rho_particle_g_cm3'), fld(meta,'rho_bulk_g_cm3'), ...
        fld(meta,'phi'), meta.fps_true, fld(meta,'nFrames'), ...
        meta.firstValidFrame, fld(meta,'nTracked'), nM, ...
        fld(meta,'windowStart'), fld(meta,'windowEnd'), ...
        double(isfield(meta,'autoWindow') && isequal(meta.autoWindow,true)));
    fclose(fid);
    fprintf('  Saved scalars CSV: %s\n', outCsv);
end

function v = fld(s, f)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end
