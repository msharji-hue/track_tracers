function D = diagnose_tracks(tracksMatPath, varargin)
% DIAGNOSE_TRACKS  Comprehensive QA on a saved _tracks.mat (read-only).
%
%   Answers, with numbers rather than impressions:
%     (Q1) Are missing markers NaN, or held at stale positions?
%     (Q2) How long is the real event window vs. the tracked window?
%     (Q3) Do markers get RE-ACQUIRED after long gaps (identity-swap risk)?
%     (Q4) Do tracked markers cross / swap order (hard ID-swap evidence)?
%     (Q5) Are inter-marker spacings physically plausible (rigid rod)?
%     (Q6) Do the "lost" markers disappear below the bed line (buried),
%          or vanish mid-air (detection failure)?
%
%   D = diagnose_tracks('/…/25mm_T01_as_poured_tracks.mat')
%   D = diagnose_tracks(path, 'detMatPath', '/…/_detections.mat', 'plot', true)
%
%   Purely diagnostic: reads files, writes an optional figure + CSV. Changes
%   nothing in the pipeline.

    p = inputParser;
    addParameter(p,'detMatPath','');
    addParameter(p,'plot',true);
    addParameter(p,'reacqGapFrames',10);   % gap beyond which re-acquisition is suspect
    parse(p,varargin{:});
    o = p.Results;

    S = load(tracksMatPath);
    X = S.tracks.trackedX;  Y = S.tracks.trackedY;
    meta = S.meta;  calib = S.calib;
    fps  = S.tracks.fps;
    [nM, nF] = size(X);
    fvf = S.tracks.firstValidFrame;

    fprintf('\n===== TRACK DIAGNOSTIC: %s =====\n', meta.trialTag);
    fprintf('markers=%d  trackedFrames=%d  fps=%.1f  firstValidFrame=%d\n', nM, nF, fps, fvf);

    % ── Q1: NaN vs stale ─────────────────────────────────────────────────
    nanMask   = ~isfinite(X);
    nanTotal  = sum(nanMask(:));
    nanPerMk  = sum(nanMask,2);
    % "stale" = identical position repeated across many consecutive frames
    stale = 0;
    for m = 1:nM
        x = X(m,:); f = isfinite(x);
        if sum(f) > 2
            d = abs(diff(x(f)));
            stale = stale + sum(d == 0);
        end
    end
    fprintf('\n[Q1] Missing-marker representation\n');
    fprintf('  NaN entries       : %d / %d (%.1f%%)\n', nanTotal, numel(X), 100*nanTotal/numel(X));
    fprintf('  exact-repeat steps: %d  (0 => no stale carry-forward)\n', stale);
    for m = 1:nM
        fprintf('    marker %d: %5d NaN (%.1f%%)\n', m, nanPerMk(m), 100*nanPerMk(m)/nF);
    end

    % ── Q2: event window vs tracked window ───────────────────────────────
    nVis   = sum(isfinite(X),1);            % visible markers per tracking frame
    full8  = find(nVis == nM);
    last8  = tern(isempty(full8), NaN, full8(end));
    % motion: frame-to-frame displacement of the leading (toe) marker = row 1
    dtoe = [NaN, sqrt(diff(X(1,:)).^2 + diff(Y(1,:)).^2)];
    movingThresh = 0.5;                     % px/frame
    moving = find(dtoe > movingThresh);
    lastMove = tern(isempty(moving), NaN, moving(end));
    fprintf('\n[Q2] Event window vs tracked window\n');
    fprintf('  frames with all %d visible : %d  (last at tracking frame %s)\n', ...
        nM, numel(full8), num2str(last8));
    fprintf('  last frame with toe motion > %.1f px : %s  (%.1f ms after track start)\n', ...
        movingThresh, num2str(lastMove), 1000*tern(isnan(lastMove),NaN,lastMove)/fps);
    fprintf('  tracked window            : %d frames (%.0f ms)\n', nF, 1000*nF/fps);
    if ~isnan(lastMove)
        fprintf('  => USEFUL fraction        : %.2f%% of tracked frames\n', 100*lastMove/nF);
    end

    % ── Q3: re-acquisition after long gaps (swap risk) ──────────────────
    fprintf('\n[Q3] Re-acquisition after gaps (>%d frames)\n', o.reacqGapFrames);
    reacq = zeros(nM,1); maxGap = zeros(nM,1);
    for m = 1:nM
        f = isfinite(X(m,:));
        d = diff(find(f));
        if ~isempty(d)
            maxGap(m) = max(d) - 1;
            reacq(m)  = sum(d-1 > o.reacqGapFrames);
        end
    end
    for m = 1:nM
        fprintf('    marker %d: max gap %5d fr (%.0f ms), long re-acquisitions: %d\n', ...
            m, maxGap(m), 1000*maxGap(m)/fps, reacq(m));
    end
    totReacq = sum(reacq);
    if totReacq > 0
        fprintf('  WARNING: %d long-gap re-acquisitions. track_markers matches against a\n', totReacq);
        fprintf('           marker''s LAST KNOWN position (any age) within tolerancePx=100.\n');
        fprintf('           Marker spacing is ~34 px, so a stale match can grab a NEIGHBOUR.\n');
    else
        fprintf('  OK: no long-gap re-acquisitions.\n');
    end

    % ── Q4: ordering / crossing (hard ID-swap evidence) ─────────────────
    % rod is rigid: x-order of markers must never change
    fprintf('\n[Q4] Marker ordering (rigid rod => x-order must be constant)\n');
    viol = 0; violFrames = [];
    for f = 1:nF
        x = X(:,f); ok = isfinite(x);
        if sum(ok) >= 2
            xs = x(ok);
            if any(diff(xs) > 0)       % row 1 = largest x (toe) => should be descending
                viol = viol + 1; violFrames(end+1) = f; %#ok<AGROW>
            end
        end
    end
    fprintf('  order-violating frames: %d / %d (%.2f%%)\n', viol, nF, 100*viol/nF);
    if viol > 0
        fprintf('  first few at tracking frames: %s\n', mat2str(violFrames(1:min(8,end))));
        fprintf('  => ID SWAP confirmed (markers crossed in x).\n');
    else
        fprintf('  OK: x-order preserved in every frame.\n');
    end

    % ── Q5: inter-marker spacing (rigid-rod sanity) ─────────────────────
    fprintf('\n[Q5] Inter-marker spacing (should be ~constant on a rigid rod)\n');
    sp = nan(nM-1, nF);
    for f = 1:nF
        x = X(:,f); y = Y(:,f);
        for m = 1:nM-1
            if isfinite(x(m)) && isfinite(x(m+1))
                sp(m,f) = hypot(x(m)-x(m+1), y(m)-y(m+1));
            end
        end
    end
    spAll = sp(isfinite(sp));
    fprintf('  spacing: mean %.2f px, sd %.2f px, min %.2f, max %.2f\n', ...
        mean(spAll), std(spAll), min(spAll), max(spAll));
    nomSp = median(spAll);
    badSp = sum(spAll > 1.6*nomSp);         % ~a skipped marker
    fprintf('  spacings > 1.6x median (skipped-marker signature): %d (%.2f%%)\n', ...
        badSp, 100*badSp/numel(spAll));

    % ── Q6: do lost markers go BELOW the bed line? ──────────────────────
    % signed perpendicular distance to bed line; sign convention as in pipeline
    nrm  = hypot(calib.lineA, calib.lineB);
    dist = (calib.lineA*X + calib.lineB*Y + calib.lineC) / nrm;
    fprintf('\n[Q6] Where do markers vanish? (bed-line signed distance)\n');
    for m = 1:nM
        f = find(isfinite(X(m,:)));
        if isempty(f), continue; end
        lastSeen = f(end);
        if lastSeen < nF
            fprintf('    marker %d last seen at frame %5d, bed-dist %7.1f px', ...
                m, lastSeen, dist(m,lastSeen));
            if dist(m,lastSeen) > 0
                fprintf('  (bed side => plausibly BURIED)\n');
            else
                fprintf('  (above bed => detection loss, NOT burial)\n');
            end
        else
            fprintf('    marker %d visible to the end\n', m);
        end
    end

    % ── pack + optional figure ───────────────────────────────────────────
    D = struct('trialTag',meta.trialTag,'nMarkers',nM,'nTracked',nF,'fps',fps, ...
        'firstValidFrame',fvf,'nanFrac',nanTotal/numel(X),'staleSteps',stale, ...
        'nFull8',numel(full8),'last8',last8,'lastMoveFrame',lastMove, ...
        'usefulFrac',tern(isnan(lastMove),NaN,lastMove/nF), ...
        'maxGap',maxGap,'reacq',reacq,'orderViolations',viol, ...
        'spacingMean',mean(spAll),'spacingSD',std(spAll),'badSpacing',badSp);

    if o.plot
        f = figure('Name','Track diagnostic','Position',[60 60 1300 800]);
        subplot(2,3,1); plot(nVis,'.-'); yline(nM,'k--'); grid on
        xlabel('tracking frame'); ylabel('# visible'); title('markers visible per frame');
        subplot(2,3,2); plot(dtoe,'.-'); yline(movingThresh,'r--'); grid on
        xlabel('tracking frame'); ylabel('toe |\Deltar| (px/frame)');
        set(gca,'YScale','log'); title('toe motion (log) — event vs rest');
        subplot(2,3,3); imagesc(isfinite(X)); colormap(gca,[1 1 1;0.2 0.5 0.9]);
        xlabel('tracking frame'); ylabel('marker'); title('visibility map (blue = tracked)');
        subplot(2,3,4); plot(sp.'); grid on
        xlabel('tracking frame'); ylabel('spacing (px)'); title('inter-marker spacing');
        subplot(2,3,5); plot(X.'); grid on
        xlabel('tracking frame'); ylabel('x (px)'); title('marker x (crossing = ID swap)');
        subplot(2,3,6); plot(dist.'); yline(0,'k--','bed'); grid on
        xlabel('tracking frame'); ylabel('signed bed distance (px)'); title('burial check');
        outPng = strrep(tracksMatPath,'.mat','_diagnostic.png');
        exportgraphics(f, outPng, 'Resolution',150);
        fprintf('\nSaved figure: %s\n', outPng);
    end

    % verdict
    fprintf('\n===== VERDICT =====\n');
    if stale == 0, fprintf('  [OK]   missing markers are NaN, not stale carry-forward\n');
    else,          fprintf('  [WARN] %d exact-repeat steps: possible stale positions\n', stale); end
    if viol == 0,  fprintf('  [OK]   no marker-order violations (no hard ID swap)\n');
    else,          fprintf('  [FAIL] %d order violations => ID SWAP present\n', viol); end
    if totReacq==0,fprintf('  [OK]   no long-gap re-acquisitions\n');
    else,          fprintf('  [WARN] %d long-gap re-acquisitions => swap risk\n', totReacq); end
    if ~isnan(lastMove) && lastMove/nF < 0.1
        fprintf('  [NOTE] only %.1f%% of tracked frames contain motion;\n', 100*lastMove/nF);
        fprintf('         the rest is the rod at rest. Consider truncating the window.\n');
    end
end

function y = tern(c,a,b), if c, y=a; else, y=b; end, end
