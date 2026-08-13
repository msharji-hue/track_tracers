function R = diag_bed(tracksPath, varargin)
% DIAG_BED  Evidence for the impact trigger: signed distance of the TOP /
% FURTHEST marker centre to the bed line, over frames, zoomed to the impact
% region — so you can see whether impactDistPx is actually reached.
%
%   R = diag_bed(tracksPath)                 % default bed line (x = 4)
%   R = diag_bed(tracksPath,'bedX',10)       % local override (trigger follows)
%   R = diag_bed(tracksPath,'zoom',300)      % frames shown around the action
%
% Prints: every marker's distance at the first full-detection frame, which one
% is furthest (the reference), that marker's full range over the trial, whether
% it CROSSES impactDistPx and at which frame, and the resulting depth.

    p = inputParser;
    addParameter(p,'bedX',[],@(x)isempty(x)||isnumeric(x));
    addParameter(p,'zoom',300,@isnumeric);
    parse(p,varargin{:}); o = p.Results;

    here = fileparts(mfilename('fullpath'));
    addpath(here, fullfile(fileparts(here),'src'));
    if isempty(o.bedX), calib = get_calibration(); else, calib = get_calibration(o.bedX); end

    S = load(tracksPath,'meta','tracks'); m = S.meta; tr = S.tracks;
    fps = m.fps_true; if ~isfinite(fps) && isfield(tr,'fps'), fps = tr.fps; end
    X = tr.trackedX; Y = tr.trackedY;

    % ── signed distances, ALL markers, ALL frames ─────────────────────────
    bedNorm = sqrt(calib.lineA^2 + calib.lineB^2);
    dAll = (calib.lineA.*X + calib.lineB.*Y + calib.lineC)./bedNorm;

    ff = find(all(isfinite(X),1),1,'first');
    if isempty(ff), ff = find(any(isfinite(X),1),1,'first'); end
    col = dAll(:,ff);
    colAbs = abs(col); colAbs(~isfinite(colAbs)) = -Inf;
    [~, refID] = max(colAbs);                    % furthest = top of rod
    ref = dAll(refID,:);
    refFin = ref(isfinite(ref));

    fprintf('\n===== %s =====\n', m.trialTag);
    fprintf('bed line x = %g px   impactDistPx = %g px   mmPerPx = %.4f\n', ...
        calib.bedX, calib.impactDistPx, calib.mmPerPx);
    fprintf('first frame with all markers detected: %d\n\n', ff);
    fprintf('marker distances at that frame (px, signed to bed line):\n');
    [~,ordr] = sort(col);
    for i = ordr(:)'
        mark = '';
        if i==refID, mark = '   <-- FURTHEST = reference (top of rod)'; end
        fprintf('   marker %d : %9.1f%s\n', i, col(i), mark);
    end
    if numel(col)>1
        sp = abs(diff(sort(col)));
        fprintf('   inter-marker spacing: median %.1f px, span %.1f px\n', ...
            median(sp,'omitnan'), max(col)-min(col));
    end

    fprintf('\nreference marker (%d) over the whole trial:\n', refID);
    fprintf('   range = [%.1f , %.1f] px\n', min(refFin), max(refFin));
    reaches = calib.impactDistPx>=min(refFin) && calib.impactDistPx<=max(refFin);
    if reaches
        fprintf('   impactDistPx = %g IS within range -> trigger is reachable\n', calib.impactDistPx);
        [~,kx] = min(abs(ref - calib.impactDistPx));
        fprintf('   closest approach at frame %d (d = %.1f px)\n', kx, ref(kx));
        sgn = sign(ref - calib.impactDistPx);
        cross = find(sgn(1:end-1).*sgn(2:end) < 0, 1, 'first');
        if ~isempty(cross)
            fprintf('   TRUE CROSSING between frames %d and %d\n', cross, cross+1);
        else
            fprintf('   (touches but does not cross)\n');
        end
    else
        fprintf('   impactDistPx = %g is OUTSIDE [%.1f, %.1f] -> not reached\n', ...
            calib.impactDistPx, min(refFin), max(refFin));
        fprintf('   nearest value the marker attains: %.1f px\n', ...
            ref(find(abs(ref-calib.impactDistPx)==min(abs(refFin-calib.impactDistPx)),1)));
    end

    % ── kinematics for context ────────────────────────────────────────────
    ws = warning('off','kd_kinematics:triggerOutOfRange');
    restore = onCleanup(@() warning(ws)); %#ok<NASGU>
    kin = kd_kinematics(X, Y, calib, 1/fps);
    fprintf('\nresulting: impact=%d  stop=%d  preFrames=%d  v0=%.1f cm/s  d=%.3f cm\n', ...
        kin.impact_index, kin.stopFrame, kin.impact_index-1, kin.v0_cm_s, kin.d_final_cm);
    % Physical height, not the label: GB/shallow labels are reversed on disk,
    % and using the label here made the free-fall ratio read 3.9-4.2.
    hTrue_mm = true_drop_height(m.dropHeight_mm, m.condition);
    v0ff     = sqrt(2*calib.g_cm_s2*(hTrue_mm/10));
    fprintf('           v0 free-fall(%g mm true) = %.1f cm/s   ratio = %.2f\n', ...
        hTrue_mm, v0ff, kin.v0_cm_s/v0ff);

    R = struct('trialTag',m.trialTag,'refID',refID,'ff',ff,'ref',ref, ...
        'range',[min(refFin) max(refFin)],'reaches',reaches,'kin',kin,'calib',calib);

    % ── plots: ZOOMED to the action, plus full-trial context ──────────────
    lo = max(1, kin.impact_index - o.zoom);
    hi = min(numel(ref), kin.stopFrame + o.zoom);

    figure('Color','w','Name',['bed check — ' m.trialTag]);

    subplot(3,1,1);
    plot(lo:hi, ref(lo:hi), '-','LineWidth',1.2); grid on; hold on;
    yline(calib.impactDistPx,'--r','impactDistPx','LabelHorizontalAlignment','left');
    xline(kin.impact_index,'-k','impact'); xline(kin.stopFrame,'-','stop');
    ylabel('top marker \rightarrow bed (px)'); xlabel('frame');
    title(sprintf('%s   ZOOMED  (marker %d = furthest)', m.trialTag, refID),'Interpreter','none');

    subplot(3,1,2);
    plot(ref,'-'); grid on; hold on;
    yline(calib.impactDistPx,'--r');
    xline(kin.impact_index,'-k'); xline(kin.stopFrame,'-');
    ylabel('same, FULL trial (px)'); xlabel('frame');
    title('full-frame view — the fall is compressed at the left');

    subplot(3,1,3);
    plot(kin.t_s*1e3, kin.z, '-'); grid on; hold on;
    xline(0,'--k','impact'); xline(kin.t_stop_s*1e3,'--r','stop');
    xlabel('t (ms, zeroed at impact)'); ylabel('depth z (cm)');
    title(sprintf('v0=%.1f cm/s   d=%.3f cm', kin.v0_cm_s, kin.d_final_cm));
end