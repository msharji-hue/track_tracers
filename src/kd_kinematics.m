function kin = kd_kinematics(trackedX, trackedY, calib, dt, varargin)
% KD_KINEMATICS  Velocity-first rod kinematics matching Katsuragi & Durian (2007).
%
%   This replaces the position-first, double-differentiated Savitzky-Golay chain
%   (smooth z -> gradient -> smooth v -> gradient -> smooth a). The literature
%   (KD 2007 Methods + Supplementary Fig. 1; corroborated structurally by the PDV
%   ballistics work) is velocity-first:
%
%       1. velocity is the primary estimate (frame-to-frame displacement);
%       2. acceleration comes from a SINGLE differentiation of velocity, via
%          ORDER-1 local straight-line fits whose window is grown until the
%          slope uncertainty falls below max(0.5%*|a|, 0.005*g)  (KD's own
%          criterion);
%       3. position is the measured rod displacement (co-registered per frame).
%
%   Our markers give position, so v requires one finite difference; a then comes
%   from ONE order-1 adaptive fit of that v. We never differentiate position
%   twice and never smooth across the stop (the acceleration is discontinuous
%   there, so segment windows are clamped to the [impact, stop] interval).
%
%   Windows are specified in MILLISECONDS and converted per trial, so the
%   physical smoothing timescale is identical across the ~2700-2800 fps spread.
%
%   INPUTS
%     trackedX/Y  [nMarkers x nFrames] tracked positions (NaN where unmatched)
%     calib       struct from get_calibration(): needs lineA/B/C, mmPerPx,
%                 impactDistPx, g_cm_s2
%     dt          1/fps_true for THIS trial (seconds)  [read per trial]
%
%   NAME-VALUE (defaults follow KD where a value is specified)
%     'minWindowMs'   min half-window for the adaptive accel fit  (default 0.5)
%     'maxWindowMs'   max half-window                              (default 4.0)
%     'accelTolRel'   relative accel-uncertainty target           (default 0.005)  % 0.5%
%     'accelTolAbs'   absolute floor, in units of g               (default 0.005)  % 0.005 g
%     'rebFactor'     stop uses v above rebFactor*max rebound      (default 2)
%     'preFallMs'     window (ms) before impact for g_eff parabola fit (default 40;
%                     needs >=~20 ms of fall for a stable per-trial g_eff)
%     'postCapMs'     ms after stop retained in z/v                (default calib.postCapMs or 10)
%
%   OUTPUT (struct kin)
%     .t_s .depthRod_cm .z .v .a .a_plus_g        per-frame series (a masked
%                                                 outside [impact,stop])
%     .impact_index .stopFrame
%     .v0_cm_s .d_final_cm .t_stop_s
%     .g_eff_cm_s2 .friction_over_m_cm_s2         pre-impact fall accel & rail friction
%     .g_eff_se_cm_s2
%     .accel_se .winHalfFrames                    per-point accel SE & window used
%     .refMarkerID .rodBedDist_px                 impact-reference metadata
%     .method                                     tag string for provenance

    p = inputParser;
    addParameter(p,'minWindowMs', 0.5,  @(x)isnumeric(x)&&x>0);
    addParameter(p,'maxWindowMs', 4.0,  @(x)isnumeric(x)&&x>0);
    addParameter(p,'accelTolRel', 0.005,@isnumeric);
    addParameter(p,'accelTolAbs', 0.005,@isnumeric);
    addParameter(p,'rebFactor',   2,    @isnumeric);
    addParameter(p,'preFallMs',   40,   @isnumeric);   % long window: g_eff scatter ~1-2 cm/s^2
    addParameter(p,'postCapMs',   [],   @(x)isempty(x)||isnumeric(x));
    parse(p,varargin{:});
    o = p.Results;

    g = calib.g_cm_s2;
    if isempty(o.postCapMs)
        if isfield(calib,'postCapMs'), o.postCapMs = calib.postCapMs; else, o.postCapMs = 10; end
    end

    % ── ms -> frames (per trial; handles the fps spread) ──────────────────
    wMin = max(1, round(o.minWindowMs*1e-3 / dt));   % half-window, frames
    wMax = max(wMin, round(o.maxWindowMs*1e-3 / dt));
    preFallFrames = max(3, round(o.preFallMs*1e-3 / dt));
    postCapFrames = round(o.postCapMs*1e-3 / dt);

    nF = size(trackedX,2);

    % ── 1) PRIMARY POSITION: composition-invariant rod displacement (cm) ──
    %     positive = into the bed. See rod_displacement.m.
    z_rod = rod_displacement(trackedX, trackedY, ...
                             calib.lineA, calib.lineB, calib.lineC, calib.mmPerPx);
    z_rod = z_rod(:).';                          % row

    % reference-marker metadata (nearest the impact plane) — not "the toe"
    bedNorm  = sqrt(calib.lineA^2 + calib.lineB^2);
    d_px_all = (calib.lineA.*trackedX + calib.lineB.*trackedY + calib.lineC)./bedNorm;
    ff = find(all(isfinite(trackedX),1),1,'first');
    if isempty(ff), ff = find(any(isfinite(trackedX),1),1,'first'); end
    col = d_px_all(:,ff); col(~isfinite(col)) = Inf;
    [~, refMarkerID] = min(abs(col - calib.impactDistPx));
    rodBedDist_px    = d_px_all(refMarkerID,:);

    t = (0:nF-1).*dt;

    % ── 2) RAW VELOCITY (one finite difference = KD's "raw v") ────────────
    v_raw = local_central_diff(z_rod, dt);

    % ── 3) IMPACT & STOP from raw velocity (no acceleration needed yet) ────
    %   Impact = peak velocity (free-fall accelerates v up to v0, drag then
    %   decelerates). Constrain the search near the geometric bed crossing so a
    %   spurious pre-release maximum can't win.
    [~, geomImpact] = min(abs(rodBedDist_px - calib.impactDistPx));
    vLite = movmean_omitnan(v_raw, 2*wMin+1);
    srch  = max(1,geomImpact-round(0.5*wMax)) : min(nF,geomImpact+round(2*wMax));
    [~, li] = max(vLite(srch));
    impact_index = srch(1) + li - 1;

    %   Stop = first zero-crossing of v after impact, refined by KD-style linear
    %   extrapolation of the segment just before the crossing, using only v above
    %   rebFactor * max post-crossing rebound speed.
    [stopFrame, t_stop_s] = find_stop(vLite, t, impact_index, o.rebFactor);

    % ── 4) ACCELERATION: order-1 adaptive line-segment fits of v_raw ──────
    %   Windows clamped to [impact, stop] so no fit straddles the stop
    %   discontinuity. v_smooth = fitted value at the point (the segment-average
    %   velocity KD plot as open circles); a = the segment slope.
    a        = nan(1,nF);
    v_smooth = nan(1,nF);
    accel_se = nan(1,nF);
    winHalf  = nan(1,nF);
    for i = 1:nF
        loEdge = 1; hiEdge = nF;
        if i >= impact_index, loEdge = impact_index; end   % don't reach pre-impact
        if i <= stopFrame,    hiEdge = stopFrame;    end    % don't reach past stop
        [a(i), v_smooth(i), accel_se(i), winHalf(i)] = ...
            adaptive_slope(t, v_raw, i, wMin, wMax, loEdge, hiEdge, ...
                           o.accelTolRel, o.accelTolAbs, g);
    end

    % net (resistive) acceleration in the project convention
    a_plus_g = -a - g;

    % ── 5) PRE-IMPACT EFFECTIVE GRAVITY (rail-friction / calibration check) ─
    %   During the guided fall the only forces are gravity and rail friction, so
    %   the free-fall obeys  z(t) = c0 + c1 t + (1/2) g_eff t^2  and g_eff = 2*c2.
    %   We fit the RAW POSITION to a parabola (not a differentiated velocity
    %   slope): position noise is iid so the quadratic coefficient is unbiased,
    %   whereas a slope-of-central-difference estimate is noise-biased. Validated
    %   on synthetic trajectories: over an 8 ms window g_eff scatter is
    %   +/-75-150 cm/s^2, but over >=40 ms it collapses to ~1-2 cm/s^2. So we use
    %   a LONG window (all available clean pre-impact frames up to preFallMs).
    pfStart = max(ff, impact_index - preFallFrames);
    pf = pfStart : (impact_index-2);
    pf = pf(pf>=1);
    if numel(pf)*dt*1e3 < 20
        warning('kd_kinematics:shortFall', ...
            ['Only %.1f ms of pre-impact fall available; g_eff is uncertain for ' ...
             'this trial (lean on the cross-trial mean).'], numel(pf)*dt*1e3);
    end
    [g_eff, g_eff_se] = fit_parabola_accel(t(pf), z_rod(pf));
    friction_over_m   = g - g_eff;          % cm/s^2   (>=0 expected)

    % ── 6) DEPTH (measured), zeroed at impact; time re-zeroed at impact ────
    depthRod_cm = z_rod - z_rod(impact_index);
    t_s         = t - t(impact_index);

    v0_cm_s    = v_smooth(impact_index);
    d_final_cm = depthRod_cm(stopFrame);

    % ── 7) MASK ────────────────────────────────────────────────────────────
    a(1:impact_index-1) = nan;               % acceleration meaningless off-interval
    a(stopFrame+1:end)  = nan;
    a_plus_g(1:impact_index-1) = nan;
    a_plus_g(stopFrame+1:end)  = nan;
    vidEnd = min(stopFrame + postCapFrames, nF);
    v_smooth(vidEnd+1:end)    = nan;
    depthRod_cm(vidEnd+1:end) = nan;

    % ── pack ────────────────────────────────────────────────────────────────
    kin = struct();
    kin.t_s           = t_s(:);
    kin.depthRod_cm   = depthRod_cm(:);
    kin.z             = depthRod_cm(:);      % alias
    kin.v             = v_smooth(:);
    kin.a             = a(:);
    kin.a_plus_g      = a_plus_g(:);
    kin.impact_index  = impact_index;
    kin.stopFrame     = stopFrame;
    kin.v0_cm_s       = v0_cm_s;
    kin.d_final_cm    = d_final_cm;
    kin.t_stop_s      = t_stop_s;
    kin.g_eff_cm_s2       = g_eff;
    kin.g_eff_se_cm_s2    = g_eff_se;
    kin.friction_over_m_cm_s2 = friction_over_m;
    kin.accel_se      = accel_se(:);
    kin.winHalfFrames = winHalf(:);
    kin.refMarkerID   = refMarkerID;
    kin.rodBedDist_px = rodBedDist_px(:);
    kin.dt            = dt;
    kin.method        = 'velocity-first; order-1 adaptive line-segment a (KD 2007)';

    fprintf(['kd_kinematics: impact=%d stop=%d | v0=%.1f cm/s d=%.3f cm t_stop=%.4f s\n' ...
             '  g_eff=%.1f cm/s^2 (%.1f%% of g) -> f/m=%.1f cm/s^2 | win %d-%d frames\n'], ...
        impact_index, stopFrame, v0_cm_s, d_final_cm, t_stop_s, ...
        g_eff, 100*g_eff/g, friction_over_m, min(winHalf(impact_index:stopFrame)), ...
        max(winHalf(impact_index:stopFrame)));
end

% ===================== helpers =========================================
function dv = local_central_diff(z, dt)
    n = numel(z); dv = nan(1,n);
    for i = 2:n-1
        if isfinite(z(i+1)) && isfinite(z(i-1))
            dv(i) = (z(i+1)-z(i-1))/(2*dt);
        elseif isfinite(z(i+1)) && isfinite(z(i))
            dv(i) = (z(i+1)-z(i))/dt;
        elseif isfinite(z(i)) && isfinite(z(i-1))
            dv(i) = (z(i)-z(i-1))/dt;
        end
    end
end

function y = movmean_omitnan(x, w)
    n = numel(x); y = nan(1,n); h = floor(w/2);
    for i = 1:n
        lo = max(1,i-h); hi = min(n,i+h);
        seg = x(lo:hi); seg = seg(isfinite(seg));
        if ~isempty(seg), y(i) = mean(seg); end
    end
end

function [slope, se] = fit_slope(x, y)
    ok = isfinite(x) & isfinite(y); x = x(ok); y = y(ok);
    n = numel(x);
    if n < 3, slope = NaN; se = NaN; return; end
    xb = mean(x); Sxx = sum((x-xb).^2);
    if Sxx == 0, slope = NaN; se = NaN; return; end
    p = polyfit(x,y,1); slope = p(1);
    resid = y - polyval(p,x);
    s2 = sum(resid.^2)/max(1,(n-2));
    se = sqrt(s2/Sxx);
end

function [accel, se] = fit_parabola_accel(x, y)
% Fit y = c2 x^2 + c1 x + c0; free-fall acceleration = 2*c2 (unbiased under iid
% position noise). SE from the coefficient covariance sigma^2 (A'A)^-1.
    ok = isfinite(x) & isfinite(y); x = x(ok); y = y(ok);
    x = x(:); y = y(:); n = numel(x);
    if n < 5, accel = NaN; se = NaN; return; end
    A = [x.^2, x, ones(n,1)];
    c = A\y;
    accel = 2*c(1);
    resid = y - A*c;
    s2 = sum(resid.^2)/max(1,(n-3));
    C  = s2 * inv(A.'*A);           %#ok<MINV>  small 3x3, fine
    se = 2*sqrt(max(0,C(1,1)));
end

function [a_i, v_i, se_i, wUsed] = adaptive_slope(t, v, i, wMin, wMax, loEdge, hiEdge, relTol, absTol, g)
% Grow a symmetric window (clamped to [loEdge,hiEdge]) until the slope SE meets
% max(relTol*|a|, absTol*g). Order-1 fit: slope=a, value-at-i=v_smooth.
    a_i = NaN; v_i = NaN; se_i = NaN; wUsed = NaN;
    for w = wMin:wMax
        lo = max(loEdge, i-w); hi = min(hiEdge, i+w);
        idx = lo:hi;
        xx = t(idx); yy = v(idx);
        ok = isfinite(xx) & isfinite(yy); xx = xx(ok); yy = yy(ok);
        if numel(xx) < 3, continue; end
        [slp, se] = fit_slope(xx, yy);
        if ~isfinite(slp), continue; end
        p = polyfit(xx,yy,1);
        a_i = slp; v_i = polyval(p, t(i)); se_i = se; wUsed = w;
        tol = max(relTol*abs(slp), absTol*g);
        if isfinite(se) && se <= tol, break; end
    end
end

function [stopFrame, t_stop] = find_stop(v, t, impact_index, rebFactor)
% First zero-crossing of v after impact, refined by linear extrapolation of the
% pre-crossing segment; ignores rebound below rebFactor*max post-crossing speed.
    n = numel(v);
    cross = [];
    for i = impact_index+1:n
        if isfinite(v(i)) && v(i) <= 0, cross = i; break; end
    end
    if isempty(cross)
        stopFrame = n; t_stop = t(n);
        warning('kd_kinematics:noStop','No v zero-crossing found; using last frame.');
        return;
    end
    reb = max(0, -min(v(cross:end),[],'omitnan'));     % rebound magnitude
    thr = rebFactor*reb;
    seg = [];
    for i = cross-1:-1:impact_index
        if isfinite(v(i)) && v(i) > thr, seg = [i seg]; else, break; end %#ok<AGROW>
    end
    if numel(seg) >= 2
        [m, ~] = fit_slope(t(seg), v(seg));
        b = mean(v(seg)) - m*mean(t(seg));
        if isfinite(m) && m ~= 0
            t_stop = -b/m;
        else
            t_stop = t(cross);
        end
    else
        t_stop = t(cross);
    end
    [~, stopFrame] = min(abs(t - t_stop));
    stopFrame = max(impact_index+1, min(stopFrame, n));
end
