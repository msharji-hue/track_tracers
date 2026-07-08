function save_tracks(meta, tracks, scalars, calib, outRoot)
% SAVE_TRACKS  Save tracked marker positions + kinematics + scalars + metadata.
%
%   Writes into  <outRoot>/tracks/ :
%       <trialTag>_tracks.mat    full bundle (meta, tracks, scalars, calib)
%       <trialTag>_tracks.csv    per-frame x/y for every marker (wide format)
%       <trialTag>_scalars.csv   one-row summary
%
%   Indexing note: trackedX/trackedY columns use TRACKING-frame numbering,
%   where column 1 == firstValidFrame in the original video. The CSV carries
%   both 'tracking_frame' and 'original_frame', plus an 'in_window' flag that
%   is 1 for frames within [1 .. savedThroughFrame] (impact + t_stop buffer).
%   The FULL (untrimmed) sequence is saved; nothing is discarded.

    trackDir = fullfile(outRoot, 'tracks');
    if ~exist(trackDir, 'dir'), mkdir(trackDir); end
    base = meta.trialTag;

    % ── Full bundle (.mat) ────────────────────────────────────────────────
    outMat = fullfile(trackDir, [base '_tracks.mat']);
    save(outMat, 'meta', 'tracks', 'scalars', 'calib', '-v7.3');
    fprintf('  Saved tracks MAT : %s\n', outMat);

    % ── Per-frame positions (.csv, wide) ─────────────────────────────────
    outCsvT = fullfile(trackDir, [base '_tracks.csv']);
    nM = size(tracks.trackedX, 1);
    nF = size(tracks.trackedX, 2);
    fid = fopen(outCsvT, 'w');

    hdr = 'tracking_frame,original_frame,t_s,in_window';
    for mm = 1:nM, hdr = [hdr sprintf(',x%d', mm)]; end %#ok<AGROW>
    for mm = 1:nM, hdr = [hdr sprintf(',y%d', mm)]; end %#ok<AGROW>
    fprintf(fid, '%s\n', hdr);

    for f = 1:nF
        origF = meta.firstValidFrame + (f - 1);
        inW   = double(f <= tracks.savedThroughFrame);
        row   = sprintf('%d,%d,%.6f,%d', f, origF, tracks.t_s(f), inW);
        for mm = 1:nM, row = [row sprintf(',%.4f', tracks.trackedX(mm,f))]; end %#ok<AGROW>
        for mm = 1:nM, row = [row sprintf(',%.4f', tracks.trackedY(mm,f))]; end %#ok<AGROW>
        fprintf(fid, '%s\n', row);
    end
    fclose(fid);
    fprintf('  Saved tracks CSV : %s\n', outCsvT);

    % ── One-row scalar summary (.csv) ────────────────────────────────────
    outCsv = fullfile(trackDir, [base '_scalars.csv']);
    fid = fopen(outCsv, 'w');
    hasBend = isfield(scalars,'bend_peak_rms_mm');
    hdr = ['material,batch,dropHeight_mm,trialNum,container,', ...
           'rho_particle_g_cm3,rho_bulk_g_cm3,phi,', ...
           'firstValidFrame,fps,v0_cm_s,d_final_cm,t_stop_s'];
    if hasBend
        hdr = [hdr, ',bend_peak_rms_mm,bend_peak_max_mm,bend_signed_pk_mm,', ...
                    'bend_at_stop_mm,bend_curv_pk_1pmm,tilt_peak_deg,', ...
                    'bend_angle_peak_deg,seg_angle_peak_deg,', ...
                    'bend_t_peak_ms,bend_baseline_mm,bend_flag,tilt_flag'];
    end
    fprintf(fid, '%s\n', hdr);
    fprintf(fid, '%s,%s,%g,%d,%s,%.4f,%.4f,%.4f,%d,%.4f,%.4f,%.4f,%.6f', ...
        meta.material, meta.batchName, meta.dropHeight_mm, meta.trialNum, ...
        meta.container, fld(meta,'rho_particle_g_cm3'), fld(meta,'rho_bulk_g_cm3'), ...
        fld(meta,'phi'), meta.firstValidFrame, meta.fps_true, ...
        scalars.v0_cm_s, scalars.d_final_cm, scalars.t_stop_s);
    if hasBend
        fprintf(fid, ',%.4f,%.4f,%.4f,%.4f,%.5g,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d', ...
            scalars.bend_peak_rms_mm, scalars.bend_peak_max_mm, ...
            scalars.bend_signed_pk_mm, scalars.bend_at_stop_mm, ...
            scalars.bend_curv_pk_1pmm, scalars.bend_tilt_pk_deg, ...
            scalars.bend_angle_pk_deg, scalars.bend_seg_pk_deg, ...
            scalars.bend_t_peak_ms, scalars.bend_baseline_mm, ...
            scalars.bend_flag, scalars.tilt_flag);
    end
    fprintf(fid, '\n');
    fclose(fid);
    fprintf('  Saved scalars CSV: %s\n', outCsv);
end

function v = fld(s, f)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end
