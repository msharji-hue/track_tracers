function diag_kd(tracksPath)
% DIAG_KD  Crash localiser for the kinematics path. Runs ONE tracks.mat with
% on-disk checkpoints, NO figures and NO dir scan. If MATLAB crashes, open the
% log — the LAST line is the step that was running when it died.
%
%   diag_kd('/full/path/to/<tag>_tracks.mat')
%
% Log: <userhome>\Desktop\diag_kd_log.txt  (LOCAL disk, not the data drive)

    home = getenv('HOME');
    if isempty(home), home = getenv('USERPROFILE'); end   % Windows
    if isempty(home), home = pwd; end
    dsk = fullfile(home,'Desktop');
    if ~exist(dsk,'dir'), dsk = home; end
    logf = fullfile(dsk,'diag_kd_log.txt');
    if exist(logf,'file'), delete(logf); end
    ck(logf, sprintf('=== diag_kd start %s ===', datestr(now)));
    ck(logf, ['tracksPath = ' tracksPath]);

    ck(logf, 'A: set path...');
    here = fileparts(mfilename('fullpath'));
    addpath(here, fullfile(fileparts(here),'src'));
    ck(logf, 'A: ok');

    ck(logf, 'B: get_calibration...');
    calib = get_calibration();
    ck(logf, sprintf('B: ok  mmPerPx=%.4f  g=%.1f', calib.mmPerPx, calib.g_cm_s2));

    ck(logf, 'C: exist(file)...');
    assert(exist(tracksPath,'file')>0, 'file not found: %s', tracksPath);
    ck(logf, 'C: ok');

    ck(logf, 'D: load meta only...');
    sm = load(tracksPath,'meta'); m = sm.meta;
    ck(logf, sprintf('D: ok  tag=%s  fps_true=%.3f', getf(m,'trialTag','?'), getf(m,'fps_true',NaN)));

    ck(logf, 'E: load tracks (this materialises the file from Dropbox)...');
    st = load(tracksPath,'tracks'); tr = st.tracks;
    ck(logf, sprintf('E: ok  size(trackedX)=%dx%d', size(tr.trackedX,1), size(tr.trackedX,2)));

    ck(logf, 'F: resolve fps/dt (resolve_fps: scalars-CSV first)...');
    [fps, fsrc] = resolve_fps(tracksPath, m, tr);
    ck(logf, sprintf('F: fps=%.3f  dt=%.6g  src=%s', fps, 1/fps, fsrc));
    assert(isfinite(fps), 'no plausible fps found (CSV/meta/tracks all out of band)');

    ck(logf, 'G: kd_kinematics (NO figure)...');
    kin = kd_kinematics(tr.trackedX, tr.trackedY, calib, 1/fps);
    % The pre-impact effective-gravity and rail-friction estimates were removed
    % from kd_kinematics: the pre-impact window (0-8 frames) is far too short to
    % fit a free-fall parabola against. a_stop is KD's acceleration
    % discontinuity, taken from the v(t) fit, and is what replaced them.
    ck(logf, sprintf('G: ok  v0=%.1f  d=%.3f  a_stop=%.1f', ...
        kin.v0_cm_s, kin.d_final_cm, kin.a_stop_cm_s2));

    ck(logf, '=== diag_kd DONE (no crash) ===');
    fprintf('\ndiag_kd completed with NO crash. Log: %s\n', logf);
end

function ck(logf, msg)
    fprintf('%s\n', msg);
    fid = fopen(logf,'a'); if fid>0, fprintf(fid,'%s\n',msg); fclose(fid); end  % survives a crash
end
function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end
end
