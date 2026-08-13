function D = diag_raw_clip(videoPath, varargin)
%DIAG_RAW_CLIP  Full pre-pipeline diagnostic on ONE raw .avi. READ-ONLY.
%
%   Answers, per clip, the four questions that decide whether the pipeline
%   can process it with stock parameters:
%     1. DURATION   how many frames / how long is the recording really?
%     2. CONTENT    is the event inside the window? (first/mid/last montage)
%     3. MARKERS    how many markers are visible, frame by frame?
%     4. DETECTION  does the EXACT pipeline chain (sharpen -> redness ->
%                   CLAHE -> median -> imfindcircles, stock params) find the
%                   expected count, and at what radius and spacing?
%
%   The radius answers the CALIBRATION question too: Default's mmPerPx was
%   set by 2 mm marker = 18.5 px (radius 9.3). Detected radii within a few
%   percent of that mean the magnification matches Default and mmPerPx
%   stands; a large deviation means this model needs its own mmPerPx.
%
%   USAGE
%       D = diag_raw_clip('D:\...\Tight Model\GB\full\0mm_T01.avi');
%       D = diag_raw_clip(p, 'nSample',40, 'FullTrace',false, 'Expected',8);
%
%   OPTIONS
%       'Expected'   expected marker count            (default 8)
%       'nSample'    frames for the detection pass    (default 40, evenly spaced)
%       'FullTrace'  blob-count EVERY frame           (default true; cheap)
%       'Montage'    show first/mid/last figure       (default true)
%
%   OUTPUT struct D: nFrames, fps_meta, dur_ms, blobCountHist, nDet, rDet,
%   spacing_px, plus verdict strings printed to the console.

    p = inputParser;
    addParameter(p,'Expected', 8);
    addParameter(p,'nSample', 40);
    addParameter(p,'FullTrace', true);
    addParameter(p,'Montage', true);
    parse(p, varargin{:}); o = p.Results;

    % stock pipeline detection parameters (mirror default_detect_params)
    P = struct('radiusRange',[9 33],'sensitivity',0.85,'edgeThresh',0.10, ...
               'polarity','bright','alphaG',0.5,'betaB',0.5,'medianK',3);
    REF_RADIUS = 18.5/2;      % px; Default calibration: 2 mm marker = 18.5 px

    [~, stem] = fileparts(videoPath);
    fprintf('\n=== diag_raw_clip : %s ===\n', stem);

    % ── 1) duration ───────────────────────────────────────────────────────
    vr = VideoReader(videoPath);
    D.nFrames  = vr.NumFrames;
    D.fps_meta = vr.FrameRate;
    D.dur_ms   = 1000*vr.Duration;
    fprintf('  frames %5d | metadata fps %.0f | duration %.0f ms | %dx%d px\n', ...
            D.nFrames, D.fps_meta, D.dur_ms, vr.Width, vr.Height);

    % ── 2) content montage ────────────────────────────────────────────────
    if o.Montage
        figure('Name', stem);
        idx = [1, round(D.nFrames/2), D.nFrames];
        lab = {'first','middle','last'};
        for k = 1:3
            subplot(3,1,k); imshow(read(vr,idx(k)));
            title(sprintf('%s (frame %d)', lab{k}, idx(k)));
        end
    end

    % ── 3) blob count, every frame (crude mask; metadata-independent) ─────
    if o.FullTrace
        n = zeros(D.nFrames,1);
        for k = 1:D.nFrames
            f = read(vr,k);
            n(k) = bwconncomp(f(:,:,1)>150 & f(:,:,2)<100).NumObjects;
        end
        D.blobCounts = n;
        fprintf('  blob count: %d frames at %d | min %d | max %d | off-count frames: %d\n', ...
                sum(n==o.Expected), o.Expected, min(n), max(n), sum(n~=o.Expected));
    end

    % ── 4) exact pipeline detection on a sample ───────────────────────────
    kk = round(linspace(1, D.nFrames, min(o.nSample, D.nFrames)));
    nDet = zeros(size(kk)); rDet = nan(size(kk)); sp = nan(size(kk));
    for j = 1:numel(kk)
        rgb = imsharpen(read(vr,kk(j)));                 % export_frames 'sharpen'
        R = im2double(rgb(:,:,1)); G = im2double(rgb(:,:,2)); B = im2double(rgb(:,:,3));
        A = max(min(R - P.alphaG*G - P.betaB*B, 1), 0);
        A = adapthisteq(A);  A = medfilt2(A,[P.medianK P.medianK]);
        [c,r] = imfindcircles(A, P.radiusRange, 'ObjectPolarity',P.polarity, ...
                  'Sensitivity',P.sensitivity, 'EdgeThreshold',P.edgeThresh);
        nDet(j) = size(c,1);
        if ~isempty(r), rDet(j) = median(r); end
        if size(c,1) >= 2
            c = sortrows(c,1);                           % rod runs along x
            sp(j) = median(diff(c(:,1)));
        end
    end
    D.nDet = nDet;  D.rDet = rDet;  D.spacing_px = median(sp,'omitnan');

    fprintf('  pipeline detects: min %d / median %g / max %d   (want %d)\n', ...
            min(nDet), median(nDet), max(nDet), o.Expected);
    fprintf('  detected radius : %.2f px  (Default calib expects %.2f; ratio %.3f)\n', ...
            median(rDet,'omitnan'), REF_RADIUS, median(rDet,'omitnan')/REF_RADIUS);
    fprintf('  marker spacing  : %.1f px  (= %.2f mm at mmPerPx 0.1079)\n', ...
            D.spacing_px, D.spacing_px*0.1079);

    % ── verdicts ──────────────────────────────────────────────────────────
    ok_det   = all(nDet == o.Expected);
    ok_scale = abs(median(rDet,'omitnan')/REF_RADIUS - 1) < 0.10;
    fprintf('\n  VERDICT detection  : %s\n', tern(ok_det, ...
        'PASS - stock parameters find the expected count on every sampled frame', ...
        'CHECK - counts deviate; inspect nDet vs frame before batch processing'));
    fprintf('  VERDICT calibration: %s\n\n', tern(ok_scale, ...
        'PASS - radius within 10%% of Default calibration; mmPerPx 0.1079 stands', ...
        'CHECK - radius differs from Default; this model may need its own mmPerPx'));
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
