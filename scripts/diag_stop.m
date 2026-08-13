function diag_stop(tightTag, defaultTag, root)
% DIAG_STOP  Why does stopFrame fire early on Tight/Wide but not Default?
%
%   Computes RAW rod displacement and velocity (no kd_kinematics masking) for
%   one trial of each model, marks the detected impact/stop, and prints the
%   frame-by-frame velocity from impact to the first zero-crossing so the
%   triggering frame is visible.
%
%   diag_stop('165mm_T03_full_tight','165mm_T03_full', 'D:\ME_GRANULAB\JerboaImpact')

    if nargin < 3, root = 'D:\ME_GRANULAB\JerboaImpact'; end
    IDP = containers.Map({'Tight','Default','Wide'}, {-376.001, -370.001, -409});
    base = get_calibration();

    tags = {tightTag, defaultTag};
    figure('Color','w','Position',[60 60 1100 850]);

    for q = 1:2
        tag = tags{q};
        D = dir(fullfile(root,'03_RESULTS','**',[tag '_tracks.mat']));
        if isempty(D), fprintf('NOT FOUND: %s\n', tag); continue; end
        tp = fullfile(D(1).folder, D(1).name);
        S  = load(tp,'meta','tracks'); m = S.meta;

        mdl = "Default";
        if endsWith(lower(string(tag)),"_tight"), mdl = "Tight";
        elseif endsWith(lower(string(tag)),"_wide"), mdl = "Wide"; end

        fps = resolve_fps(tp, m, S.tracks);
        dt  = 1/fps;
        cb  = base; cb.impactDistPx = IDP(char(mdl));

        % ---- RAW, unmasked ----
        zr = rod_displacement(S.tracks.trackedX, S.tracks.trackedY, ...
                              cb.lineA, cb.lineB, cb.lineC, cb.mmPerPx);
        zr = zr(:).';
        nF = numel(zr);
        vr = nan(1,nF);
        for i = 2:nF-1
            if isfinite(zr(i+1)) && isfinite(zr(i-1)), vr(i) = (zr(i+1)-zr(i-1))/(2*dt); end
        end

        % ---- what kd_kinematics decides ----
        ws = warning('off','all');
        evalc('kin = kd_kinematics(S.tracks.trackedX,S.tracks.trackedY,cb,dt);');
        warning(ws);
        i0 = kin.impact_index; sf = kin.stopFrame;

        % ---- where would z stop rising if we ignored the stop rule? ----
        zi = zr - zr(i0);
        [zmax, imax] = max(zi(i0:min(nF, i0+round(0.100/dt))));   % within 100 ms
        imax = imax + i0 - 1;

        fprintf('\n================ %s  (%s) ================\n', tag, mdl);
        fprintf('fps=%.0f  dt=%.3f ms  nFrames=%d\n', fps, dt*1e3, nF);
        fprintf('impact=%d  stop=%d   -> nPost = %d frames = %.2f ms\n', ...
                i0, sf, sf-i0, (sf-i0)*dt*1e3);
        fprintf('d_final (at stop)          = %.4f cm\n', kin.d_final_cm);
        fprintf('MAX raw z within 100 ms    = %.4f cm  at frame %d (%.2f ms after impact)\n', ...
                zmax, imax, (imax-i0)*dt*1e3);
        fprintf('  --> if the stop rule is truncating, these two differ a lot.\n');

        % ---- the frames that trigger the crossing ----
        fprintf('\nframe-by-frame raw v from impact to first v<=0:\n');
        shown = 0;
        for i = i0:min(nF, i0+40)
            fprintf('   f=%-6d t=%+7.3f ms   v=%+10.1f cm/s   z=%+7.4f cm', ...
                    i, (i-i0)*dt*1e3, vr(i), zi(i));
            if isfinite(vr(i)) && vr(i) <= 0
                fprintf('   <== FIRST v<=0 : cross fires here');
                fprintf('\n'); shown = 1; break;
            end
            fprintf('\n');
        end
        if ~shown, fprintf('   (no v<=0 within 40 frames of impact)\n'); end

        % ---- plots ----
        tms = ((1:nF)-i0)*dt*1e3;
        subplot(2,2,q); hold on; grid on;
        plot(tms, zi, '-', 'LineWidth',1);
        xline(0,'k-','impact'); xline((sf-i0)*dt*1e3,'r--','stop');
        plot((imax-i0)*dt*1e3, zmax, 'ko','MarkerFaceColor','y','MarkerSize',7);
        xlabel('t (ms, 0 = impact)'); ylabel('RAW z (cm), unmasked');
        title(sprintf('%s — %s   d@stop=%.3f  zmax=%.3f cm', mdl, tag, kin.d_final_cm, zmax), ...
              'Interpreter','none');
        xlim([-20 100]);

        subplot(2,2,q+2); hold on; grid on;
        plot(tms, vr, '-', 'LineWidth',1);
        yline(0,'k:'); xline(0,'k-'); xline((sf-i0)*dt*1e3,'r--','stop');
        xlabel('t (ms, 0 = impact)'); ylabel('RAW v (cm/s)');
        title('raw velocity (no smoothing)'); xlim([-20 100]);
    end
end
