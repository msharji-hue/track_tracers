function sg_sensitivity(rawMatPath, varargin)
% SG_SENSITIVITY  Before/after diagnostic for the SG smoothing/differentiation
% method, run on ONE trial's raw tracks (…_tracks_raw.mat). Does NOT touch the
% pipeline — read-only exploration so you can choose defaults before locking.
%
%   sg_sensitivity('/…/25mm_T01_full_tracks_raw.mat')
%   sg_sensitivity(path, 'orders',[2 3 4], 'windowsMs',[3 5 8 10], ...
%                        'defOrder',3, 'defWinMs',5, 'prerollMs',8, 'postrollMs',8)
%
%   Produces:
%     (1) a 6-panel before/after figure (OLD full-window method vs NEW default),
%     (2) a sweep table over orders x windows (v0, t_stop, d_final, negative-a+g
%         fraction, and the a+g-vs-v^2 linear-fit R^2), printed and saved to CSV.
%
%   OLD method = current toe_kinematics (full-signal SG, window = 10% of valid,
%   compounded sgolayfilt-of-gradient, stop at a+g<=0).
%   NEW method = impact-window-confined single-fit SG (sg_derivatives), stop at
%   first v<=0 after impact (max penetration).

    p = inputParser;
    addParameter(p,'orders',[2 3 4]);
    addParameter(p,'windowsMs',[3 5 8 10]);
    addParameter(p,'defOrder',3);
    addParameter(p,'defWinMs',5);
    addParameter(p,'prerollMs',8);
    addParameter(p,'postrollMs',8);
    addParameter(p,'vFitFrac',0.15);   % exclude tail below this fraction of v0
    parse(p,varargin{:});
    o = p.Results;

    % ── load raw tracks + calibration ────────────────────────────────────
    S = load(rawMatPath);
    trackedX = S.trackedX;  trackedY = S.trackedY;  fps = S.fps;
    dt = 1/fps;  g = 980;
    c  = get_calibration();

    % ── toe depth series (same projection as toe_kinematics; read-only) ──
    ff = find(any(isfinite(trackedX),1),1,'first');
    [~,toeID] = max(trackedX(:,ff),[],'omitnan');
    bedNorm = sqrt(c.lineA^2 + c.lineB^2);
    toePx = (c.lineA.*trackedX(toeID,:) + c.lineB.*trackedY(toeID,:) + c.lineC)./bedNorm;
    [~,impact] = min(abs(toePx - c.impactDistPx));
    nF = numel(toePx);
    trel = ((0:nF-1)*dt) - (impact-1)*dt;
    depth = (toePx - toePx(impact)).*c.mmPerPx./10;    % cm
    valid = isfinite(depth);
    vIdx  = find(valid);  zc = depth(vIdx).';  nValid = numel(zc);

    % ═════════ OLD method (current toe_kinematics) ═════════
    wOld = max(7, round(0.10*nValid)); if mod(wOld,2)==0, wOld=wOld+1; end
    zO = sgolayfilt(zc, 2, wOld);
    vO = sgolayfilt(gradient(zO,dt), 2, wOld);
    aO = sgolayfilt(gradient(vO,dt), 2, wOld);
    z_old=nan(nF,1); v_old=nan(nF,1); a_old=nan(nF,1);
    z_old(vIdx)=zO; v_old(vIdx)=vO; a_old(vIdx)=aO;
    ag_old = -a_old - g;
    so = find(ag_old(impact:end)<=0,1,'first'); stop_old = tern(isempty(so),nF,so+impact-1);

    % ═════════ NEW method: provisional stop to bound the window ═════════
    [~, vprov, ~] = sg_derivatives(zc, dt, 2, oddf(round(0.008*fps)));
    vprovF = nan(nF,1); vprovF(vIdx)=vprov;
    sp = find(vprovF(impact:end)<=0,1,'first'); stopProv = tern(isempty(sp),nF,sp+impact-1);
    preF  = round(o.prerollMs/1000*fps);
    postF = round(o.postrollMs/1000*fps);
    w0 = max(1, impact-preF);  w1 = min(nF, stopProv+postF);
    winIdx = (w0:w1).';                         % impact window (frames)

    % ═════════ sweep table ═════════
    rows = {};
    fprintf('\n order  win_ms  win_fr   v0(cm/s)  t_stop(ms)  d_fin(cm)  negAG%%  R2(a+g~v2)\n');
    for ord = o.orders
        for wm = o.windowsMs
            R = new_kinematics(depth, impact, dt, fps, ord, wm, winIdx, g, o.vFitFrac);
            rows(end+1,:) = {ord, wm, R.win, R.v0, R.tstop_ms, R.dfin, R.negfrac*100, R.r2}; %#ok<AGROW>
            fprintf('  %d      %2d     %3d     %7.1f    %7.2f     %6.3f   %5.1f    %.4f\n', ...
                ord, wm, R.win, R.v0, R.tstop_ms, R.dfin, R.negfrac*100, R.r2);
        end
    end
    outCsv = strrep(rawMatPath,'.mat','_sgsweep.csv');
    write_sweep(outCsv, rows);
    fprintf('Saved sweep: %s\n', outCsv);

    % ═════════ before/after figure (OLD vs NEW default) ═════════
    D = new_kinematics(depth, impact, dt, fps, o.defOrder, o.defWinMs, winIdx, g, o.vFitFrac);

    f = figure('Name','SG before/after','Position',[60 60 1300 820]);
    tl = @(s) title(s,'FontWeight','bold');
    tw = tern(isfinite(stop_old), max(stop_old,D.stop)+postF, nF);
    xr = [max(1,impact-preF), min(nF,tw)];

    subplot(2,3,1); hold on
    plot(trel(vIdx)*1000, depth(vIdx), '.', 'Color',[.7 .7 .7]);
    plot(trel*1000, z_old, 'r', 'LineWidth',1.2);
    plot(trel(D.gi)*1000, D.z, 'b', 'LineWidth',1.2);
    xlim(trel(xr)*1000); grid on; xlabel('t - t_{impact} (ms)'); ylabel('z (cm)');
    legend({'raw','OLD','NEW'},'Location','best'); tl('depth z(t)');

    subplot(2,3,2); hold on
    plot(trel*1000, v_old, 'r'); plot(trel(D.gi)*1000, D.v, 'b');
    yline(0,'k:'); xlim(trel(xr)*1000); grid on
    xlabel('t - t_{impact} (ms)'); ylabel('v (cm/s)'); tl('velocity v(t)');

    subplot(2,3,3); hold on
    plot(trel*1000, ag_old, 'r'); plot(trel(D.gi)*1000, D.ag, 'b');
    yline(0,'k:');
    xline(trel(stop_old)*1000,'r--','stop_{old}');
    xline(trel(D.stop)*1000,'b--','stop_{new}');
    xlim(trel(xr)*1000); grid on
    xlabel('t - t_{impact} (ms)'); ylabel('a+g (cm/s^2)'); tl('a+g(t)  (fit region shaded)');
    yl=ylim; patch(trel([D.fitS D.fitE D.fitE D.fitS])*1000,[yl(1) yl(1) yl(2) yl(2)], ...
        [0 0 1],'FaceAlpha',0.06,'EdgeColor','none');

    subplot(2,3,4); hold on
    plot(v_old.^2, ag_old, '.', 'Color',[1 .5 .5]);
    plot(D.v.^2, D.ag, '.', 'Color',[.4 .4 1]);
    yline(0,'k:'); grid on; xlabel('v^2 (cm/s)^2'); ylabel('a+g (cm/s^2)');
    legend({'OLD','NEW'},'Location','best'); tl('a+g vs v^2 (all)');

    subplot(2,3,5); hold on
    plot(D.vfit.^2, D.agfit, 'b.', 'MarkerSize',8);
    pf = polyfit(D.vfit.^2, D.agfit, 1); xx=linspace(min(D.vfit.^2),max(D.vfit.^2),50);
    plot(xx, polyval(pf,xx), 'k-'); grid on
    xlabel('v^2 (cm/s)^2'); ylabel('a+g (cm/s^2)');
    tl(sprintf('NEW fit region: slope=%.4g  R^2=%.4f', pf(1), D.r2));

    subplot(2,3,6); axis off
    txt = {
        sprintf('OLD: win=%d fr (%.0f ms), stop=%d', wOld, wOld*dt*1000, stop_old)
        sprintf('NEW: order=%d win=%d fr (%.0f ms)', o.defOrder, D.win, D.win*dt*1000)
        sprintf('     preroll=%.0f ms  postroll=%.0f ms', o.prerollMs, o.postrollMs)
        sprintf('impact frame = %d', impact)
        sprintf('stop_old=%d   stop_new=%d', stop_old, D.stop)
        sprintf('v0(NEW)=%.1f cm/s', D.v0)
        sprintf('t_stop(NEW)=%.2f ms', D.tstop_ms)
        sprintf('d_final(NEW)=%.3f cm', D.dfin)
        sprintf('neg a+g fraction (NEW)=%.1f%%', D.negfrac*100)
        sprintf('fit v-threshold=%.0f%% of v0', o.vFitFrac*100) };
    text(0.02,0.98,txt,'VerticalAlignment','top','FontName','FixedWidth','FontSize',10,'Interpreter','none');
    tl('settings');

    outPng = strrep(rawMatPath,'.mat','_sgdiag.png');
    exportgraphics(f, outPng, 'Resolution',150);
    fprintf('Saved figure: %s\n', outPng);
end

% ── NEW-method kinematics on the impact window, returns metrics ────────────
function R = new_kinematics(depth, impact, dt, fps, ord, wm, winIdx, g, vFitFrac)
    win = oddf(round(wm/1000*fps));
    zwin = depth(winIdx).';
    gm   = isfinite(zwin);
    zc   = zwin(gm);
    gidxLocal = find(gm);
    gi   = winIdx(gm);                      % global frame indices (finite)
    [zs,vs,as] = sg_derivatives(zc, dt, ord, win);
    ag = -as - g;

    impLocal = find(gi==impact,1);
    if isempty(impLocal), impLocal = max(1, find(gi>=impact,1)); end

    % stop = first v<=0 after impact (max penetration)
    sIdx = find(vs(impLocal:end) <= 0, 1, 'first');
    if isempty(sIdx), stopLocal = numel(zc); else, stopLocal = sIdx+impLocal-1; end
    stop = gi(stopLocal);

    v0 = vs(impLocal);
    dfin = zs(stopLocal);
    tstop_ms = (gi(stopLocal)-impact)*dt*1000;

    % fit region: impact(+1) .. last frame with v >= vFitFrac*v0, capped at stop
    vth = vFitFrac*abs(v0);
    fitMask = false(size(vs));
    fitMask(impLocal:stopLocal) = abs(vs(impLocal:stopLocal)) >= vth;
    fS_local = find(fitMask,1,'first'); fE_local = find(fitMask,1,'last');
    if isempty(fS_local), fS_local=impLocal; fE_local=stopLocal; end
    vfit = vs(fS_local:fE_local); agfit = ag(fS_local:fE_local);
    ok = isfinite(vfit)&isfinite(agfit);
    vfit=vfit(ok); agfit=agfit(ok);
    if numel(vfit)>=3
        pf = polyfit(vfit.^2, agfit, 1); yp=polyval(pf,vfit.^2);
        r2 = 1 - sum((agfit-yp).^2)/max(eps,sum((agfit-mean(agfit)).^2));
    else, r2 = NaN; end

    segStop = ag(impLocal:stopLocal);
    negfrac = mean(segStop < 0, 'omitnan');

    R = struct('win',win,'v0',v0,'tstop_ms',tstop_ms,'dfin',dfin,'r2',r2, ...
               'negfrac',negfrac,'z',zs,'v',vs,'ag',ag,'gi',gi,'stop',stop, ...
               'v0_lbl',v0,'fitS',gi(fS_local),'fitE',gi(fE_local), ...
               'vfit',vfit,'agfit',agfit);
    R.v0 = v0;
end

function w = oddf(w), w=round(w); if mod(w,2)==0, w=w+1; end, if w<5, w=5; end, end
function y = tern(c,a,b), if c, y=a; else, y=b; end, end

function write_sweep(path, rows)
    fid=fopen(path,'w');
    fprintf(fid,'sg_order,sg_window_ms,sg_window_frames,v0_cm_s,t_stop_ms,d_final_cm,neg_ag_pct,r2_ag_v2\n');
    for i=1:size(rows,1)
        fprintf(fid,'%d,%d,%d,%.3f,%.3f,%.4f,%.2f,%.4f\n', rows{i,1},rows{i,2},rows{i,3}, ...
            rows{i,4},rows{i,5},rows{i,6},rows{i,7},rows{i,8});
    end
    fclose(fid);
end
