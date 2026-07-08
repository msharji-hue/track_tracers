function bend = rod_bending(trackedX, trackedY, impact_index, stopFrame, t_s, mmPerPx, varargin)
% ROD_BENDING  Quantify transient rod flexure from the tracked marker centerline.
%
%   The 8 circle centers sample the rod axis. In every frame the axis is fit by
%   TOTAL LEAST SQUARES (PCA), which is orientation-free (works for a near-
%   vertical rod, unlike OLS y=f(x)). PC1 = rod axis, the projection onto PC2 =
%   signed lateral deflection of each marker. A pre-impact BASELINE shape is
%   subtracted per marker so we measure CHANGE in shape (bending), not static
%   non-straightness. Rigid TILT (rotation of PC1 vs. reference) is reported
%   separately, since rotation is not bending.
%
%   bend = rod_bending(trackedX, trackedY, impact_index, stopFrame, t_s, mmPerPx, ...)
%
%   Inputs:
%     trackedX,trackedY : [nMarkers x nFrames] tracked positions (tracking-frame idx)
%     impact_index      : frame index of impact (into the tracked arrays)
%     stopFrame         : frame index of stop
%     t_s               : [1 x nFrames] time vector (s) for peak timing
%     mmPerPx           : calibration (mm per pixel)
%   Name-value:
%     'preWindow'   (default 8)     # pre-impact frames used for the baseline shape
%     'relFactor'   (default 3)     bendFlag if peak_rms > relFactor*baseline_rms
%     'absThresh_mm'(default 0.10)  ... and peak_rms > this floor (mm)
%     'minMarkers'  (default 3)     min valid markers for a frame to be scored
%
%   Output struct bend (per-frame time series + scalar summaries), all bending
%   lengths in mm (…_px twins kept for provenance).

    p = struct('preWindow',8, 'relFactor',3, 'absThresh_mm',0.10, 'minMarkers',3, ...
               'angleThresh_deg',1.0, 'tiltThresh_deg',2.0);
    for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end

    [nM, nF] = size(trackedX);

    % ── pass 1: raw per-frame PCA (axis, centre) ─────────────────────────
    PC1 = nan(2,nF); PC2 = nan(2,nF); Cen = nan(2,nF); nValid = zeros(1,nF);
    for f = 1:nF
        P  = [trackedX(:,f), trackedY(:,f)];
        ok = all(isfinite(P),2);
        nValid(f) = sum(ok);
        if nValid(f) < p.minMarkers, continue; end
        c  = mean(P(ok,:),1);
        Q  = P(ok,:) - c;
        [~,~,V] = svd(Q, 0);
        PC1(:,f) = V(:,1); PC2(:,f) = V(:,2); Cen(:,f) = c(:);
    end

    % ── baseline window (anchored to impact; fallback to leading frames) ──
    lowConfidence = false;
    b0 = max(1, impact_index - p.preWindow);
    b1 = max(1, impact_index - 1);
    bwin = b0:b1;
    bwin = bwin(nValid(bwin) >= p.minMarkers);
    if numel(bwin) < 2                       % too few clean pre-impact frames
        cand = find(nValid >= p.minMarkers, min(p.preWindow,nF), 'first');
        bwin = cand(:)'; lowConfidence = true;
    end

    % reference orientation = baseline frame with the most valid markers
    [~,mi] = max(nValid(bwin)); refIdx = bwin(mi);
    refPC1 = PC1(:,refIdx); refPC2 = PC2(:,refIdx);

    % ── pass 2: oriented residuals e(marker,frame) & tilt ────────────────
    E = nan(nM,nF); S = nan(nM,nF); tilt_deg = nan(1,nF);
    for f = 1:nF
        if any(isnan(PC1(:,f))), continue; end
        u = PC1(:,f); w = PC2(:,f);
        if dot(u,refPC1) < 0, u = -u; end     % consistent axis direction
        if dot(w,refPC2) < 0, w = -w; end     % consistent normal direction
        P  = [trackedX(:,f), trackedY(:,f)];
        ok = all(isfinite(P),2);
        Q  = P - Cen(:,f)';                   % centre (NaNs stay NaN)
        S(ok,f) = Q(ok,:) * u;                % along-axis coord (px)
        E(ok,f) = Q(ok,:) * w;                % signed perpendicular dev (px)
        tilt_deg(f) = acosd(min(1, abs(dot(u,refPC1))));
    end

    % ── baseline shape (per marker) and bending residual ─────────────────
    baseline_e = nmean(E(:,bwin), 2);             % static non-straightness, px
    Eb = E - baseline_e;                          % bending residual, px
    tmpB = E(:,bwin).^2;
    baseline_rms_px = sqrt(nmean(tmpB(:), 1));

    % ── per-frame metrics ────────────────────────────────────────────────
    rms_px = sqrt(nmean(Eb.^2, 1));
    absEb = abs(Eb); absEb(isnan(absEb)) = -Inf;
    [maxAbs_px, mrow] = max(absEb, [], 1);
    allnan = all(isnan(Eb), 1);
    maxAbs_px(allnan) = NaN;
    mrow = double(mrow); mrow(allnan) = NaN;
    signedMax_px   = nan(1,nF);
    curv_1ppx      = nan(1,nF);
    bend_angle_deg = nan(1,nF);   % upper-half vs lower-half segment angle
    seg_angle_deg  = nan(1,nF);   % max local neighbour-segment kink
    for f = 1:nF
        if isfinite(mrow(f)), signedMax_px(f) = Eb(mrow(f),f); end
        ok = isfinite(S(:,f)) & isfinite(Eb(:,f));
        if sum(ok) < 3, continue; end
        sO = S(ok,f); eO = Eb(ok,f);
        [sO, ordr] = sort(sO); eO = eO(ordr);    % order markers along the rod
        % curvature from a quadratic fit of deflection vs along-axis
        a = polyfit(sO, eO, 2);
        curv_1ppx(f) = 2*a(1);
        % angular metrics are computed in the rod frame (tilt-invariant,
        % baseline already removed via Eb) so a straight/tilted rod reads 0
        if numel(sO) >= 4
            med = median(sO); lo = sO <= med; hi = ~lo;
            if sum(lo) >= 2 && sum(hi) >= 2
                pl = polyfit(sO(lo), eO(lo), 1);
                ph = polyfit(sO(hi), eO(hi), 1);
                bend_angle_deg(f) = atand(ph(1)) - atand(pl(1));
            end
        end
        if numel(sO) >= 3
            segSlope = diff(eO) ./ diff(sO);
            dseg = diff(atand(segSlope));         % kink at each interior joint
            if ~isempty(dseg)
                [~,ik] = max(abs(dseg)); seg_angle_deg(f) = dseg(ik);
            end
        end
    end

    % ── convert to mm ────────────────────────────────────────────────────
    rms_mm       = rms_px       * mmPerPx;
    max_mm       = maxAbs_px    * mmPerPx;
    signedMax_mm = signedMax_px * mmPerPx;
    curv_1pmm    = curv_1ppx    / mmPerPx;        % 1/px -> 1/mm
    baseline_rms_mm = baseline_rms_px * mmPerPx;

    % ── summarise over the impact -> stop window ─────────────────────────
    w0 = max(1, impact_index);
    w1 = min(nF, max(stopFrame, w0));
    win = w0:w1;

    [peak_rms_mm, ir] = maxfin(rms_mm(win));
    peakFrame = win(ir);
    peak_max_mm    = maxfin(max_mm(win));
    signed_peak_mm = signedMax_mm(peakFrame);
    bend_at_stop_mm= valat(rms_mm, min(stopFrame,nF));
    [~, ic] = maxfin(abs(curv_1pmm(win)));
    curv_peak_1pmm = curv_1pmm(win(ic));
    tilt_peak_deg  = maxfin(tilt_deg(win));
    [~, iba] = maxfin(abs(bend_angle_deg(win)));
    if isfinite(iba), bend_angle_peak_deg = bend_angle_deg(win(iba)); else, bend_angle_peak_deg = NaN; end
    seg_angle_peak_deg = maxfin(abs(seg_angle_deg(win)));

    if isfinite(peakFrame) && isfinite(impact_index) && impact_index>=1 && ...
       peakFrame<=numel(t_s) && impact_index<=numel(t_s)
        t_peak_ms = (t_s(peakFrame) - t_s(impact_index)) * 1000;
    else
        t_peak_ms = NaN;
    end

    mmHit  = isfinite(peak_rms_mm) && ...
             peak_rms_mm > max(p.absThresh_mm, p.relFactor*baseline_rms_mm);
    angHit = isfinite(bend_angle_peak_deg) && abs(bend_angle_peak_deg) > p.angleThresh_deg;
    bendFlag = mmHit || angHit;
    tiltFlag = isfinite(tilt_peak_deg) && tilt_peak_deg > p.tiltThresh_deg;

    % ── pack ─────────────────────────────────────────────────────────────
    bend = struct();
    bend.method          = 'PCA/TLS perpendicular deflection, baseline-subtracted';
    bend.mmPerPx         = mmPerPx;
    % time series
    bend.rms_mm          = rms_mm;
    bend.max_mm          = max_mm;
    bend.signedMax_mm    = signedMax_mm;
    bend.curv_1pmm       = curv_1pmm;
    bend.tilt_deg        = tilt_deg;
    bend.bend_angle_deg  = bend_angle_deg;
    bend.seg_angle_deg   = seg_angle_deg;
    bend.residual_px     = Eb;              % [nM x nF] signed, for deep dives
    % baseline / reference
    bend.baseline_e_px   = baseline_e(:)';
    bend.baseline_rms_mm = baseline_rms_mm;
    bend.baselineFrames  = bwin;
    bend.refFrame        = refIdx;
    bend.lowConfidence   = lowConfidence;
    % summaries (scalars)
    bend.peak_rms_mm     = peak_rms_mm;
    bend.peak_max_mm     = peak_max_mm;
    bend.signed_peak_mm  = signed_peak_mm;
    bend.bend_at_stop_mm = bend_at_stop_mm;
    bend.curv_peak_1pmm  = curv_peak_1pmm;
    bend.tilt_peak_deg   = tilt_peak_deg;
    bend.bend_angle_peak_deg = bend_angle_peak_deg;
    bend.seg_angle_peak_deg  = seg_angle_peak_deg;
    bend.t_peak_ms       = t_peak_ms;
    bend.peakFrame       = peakFrame;
    bend.bendFlag        = bendFlag;
    bend.tiltFlag        = tiltFlag;
    bend.relFactor       = p.relFactor;
    bend.absThresh_mm    = p.absThresh_mm;
    bend.angleThresh_deg = p.angleThresh_deg;
    bend.tiltThresh_deg  = p.tiltThresh_deg;
end

% ── helpers ───────────────────────────────────────────────────────────────
function m = nmean(X, dim)
    n  = sum(~isnan(X), dim);
    X0 = X; X0(isnan(X0)) = 0;
    m  = sum(X0, dim) ./ n;
    m(n==0) = NaN;
end

function [v,i] = maxfin(x)
    x = x(:)'; m = isfinite(x);
    if ~any(m), v = NaN; i = NaN; return; end
    idx = find(m); [v,k] = max(x(idx)); i = idx(k);
end

function v = valat(x, i)
    if isfinite(i) && i>=1 && i<=numel(x), v = x(i); else, v = NaN; end
end
