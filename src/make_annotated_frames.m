function make_annotated_frames(framesDir, detectOut, trackedX, trackedY, ...
                               firstValidFrame, calib, outDir, varargin)
% MAKE_ANNOTATED_FRAMES  Write annotated PNGs starting at firstValidFrame.
%
%   Every output frame gets:
%     • dashed YELLOW bed line (from calib.bedPoint1/2, drawn across the frame)
%     • thin circle outline at each detection (shows detected radius/quality)
%     • thick centre dot at each TRACKED marker, coloured by marker ID
%   Per-ID colours make an ID swap immediately visible — the main tracking
%   failure mode. Uniform colour would hide it.
%
%   Frame k of the output corresponds to EXPORTED-WINDOW frame
%   (firstValidFrame + k - 1) and TRACKING frame k. The exported window is not
%   necessarily the whole video (see opts.autoWindow in process_trial), so the
%   absolute video frame is meta.windowStart + firstValidFrame + k - 2.
%
%   make_annotated_frames(framesDir, detectOut, trackedX, trackedY, ...
%                         firstValidFrame, calib, outDir, ...
%                         'scale',6, 'showIDs',false, 'maxFrames',Inf)
%
%   Inputs
%     framesDir       folder of exported frame_*.png (window indexing)
%     detectOut       struct from detect_circles_per_frame (centersCell/radiiCell)
%     trackedX/Y      [nMarkers x nTrackFrames], column 1 == firstValidFrame
%     firstValidFrame window-frame index that became tracking frame 1
%     calib           get_calibration() struct (bedPoint1/bedPoint2)
%     outDir          destination for annotated PNGs
%
%   Options
%     'scale'      integer upsample factor for legibility (default 6; the raw
%                  frame is ~416x36 px, far too small to read at 1:1)
%     'showIDs'    draw marker-ID numerals (default false; cramped at 36 px)
%     'maxFrames'  cap the number of annotated frames written (default Inf)
%     'Video'      path to the raw video. Supply this (with framesDir = '')
%                  when Stage A streamed and wrote no PNGs: the renderer then
%                  re-reads only the frames it actually needs, straight from
%                  the video, instead of requiring a 01_FRAMES cache.
%     'WindowStart' absolute video frame of window index 1. Required with
%                  'Video', since window index k is video frame
%                  WindowStart + k - 1.
%     'Filter'     filter to apply to video-read frames, matching the one Stage
%                  A used. The PNG cache stored filtered frames, so without
%                  this the overlays would differ cosmetically between the two
%                  paths.

    p = inputParser;
    addParameter(p,'scale',6);
    addParameter(p,'showIDs',false);
    addParameter(p,'maxFrames',Inf);
    addParameter(p,'Video','');
    addParameter(p,'WindowStart',1);
    addParameter(p,'Filter','');
    parse(p,varargin{:});
    S       = max(1, round(p.Results.scale));
    showIDs = p.Results.showIDs;

    if ~exist(outDir,'dir'), mkdir(outDir); end

    % Frame provider: either the PNG cache or the raw video. Both are indexed
    % by WINDOW index, so everything below is identical either way.
    useVideo = ~isempty(p.Results.Video);
    if useVideo
        vr      = VideoReader(p.Results.Video);
        wStart  = p.Results.WindowStart;
        nAvail  = floor(vr.Duration * vr.FrameRate);
        nInWin  = nAvail - wStart + 1;
        % The PNG cache held FILTERED frames, so apply the same filter here or
        % the overlays would look different depending on how Stage A ran.
        filt    = p.Results.Filter;
        readWin = @(k) local_read_window_frame(vr, wStart, k, filt);
    else
        files   = dir(fullfile(framesDir,'*.png'));
        files   = sort({files.name});
        nInWin  = numel(files);
        readWin = @(k) imread(fullfile(framesDir, files{k}));
    end

    nTrack = size(trackedX,2);
    nOut   = min([nTrack, nInWin-firstValidFrame+1, p.Results.maxFrames]);
    if nOut < 1, error('make_annotated_frames: nothing to render.'); end

    nM   = size(trackedX,1);
    cmap = marker_palette(nM);          % fixed palette; no dependence on lines()

    % bed line endpoints in ORIGINAL pixel coords, extended across the frame
    probe = readWin(firstValidFrame);
    H = size(probe,1); W = size(probe,2);
    [bx, by] = bed_endpoints(calib, W, H);

    for k = 1:nOut
        winIdx = firstValidFrame + k - 1;
        img = readWin(winIdx);
        if size(img,3)==1, img = repmat(img,1,1,3); end
        img = imresize(img, S, 'nearest');           % upsample for legibility

        % ── dashed yellow bed line ───────────────────────────────────────
        img = draw_dashed_line(img, bx*S, by*S, uint8([255 230 50]), max(2,round(S/2)));

        % ── thin outline at each DETECTION (radius quality) ──────────────
        c = detectOut.centersCell{winIdx};
        r = detectOut.radiiCell{winIdx};
        for i = 1:size(c,1)
            img = draw_circle(img, c(i,1)*S, c(i,2)*S, r(i)*S, uint8([255 255 255]), 2);
        end

        % ── thick centre dot per TRACKED marker, coloured by ID ──────────
        for m = 1:nM
            xm = trackedX(m,k); ym = trackedY(m,k);
            if ~isfinite(xm) || ~isfinite(ym), continue; end
            col = uint8(round(cmap(m,:)*255));
            img = draw_disc(img, xm*S, ym*S, max(3,round(S*1.4)), col);
            if showIDs
                img = insertText(img, [xm*S, ym*S - 4*S], sprintf('%d',m), ...
                    'FontSize',max(8,2*S), 'BoxOpacity',0, 'TextColor',col);
            end
        end

        % ── warn when the detected count is not 8 ────────────────────────
        nDet = size(c,1);
        if nDet ~= 8
            img = draw_border(img, uint8([255 40 40]), max(2,round(S)));
        end

        imwrite(img, fullfile(outDir, sprintf('annot_%05d.png', k)));
    end
    fprintf('  Annotated %d frames (orig %d..%d) -> %s\n', ...
        nOut, firstValidFrame, firstValidFrame+nOut-1, outDir);
end

% ── helpers (pure pixel drawing; no graphics/toolbox colormap deps) ────────
function cmap = marker_palette(n)
    base = [ 0.20 0.90 0.20     % green
             0.20 0.70 1.00     % cyan-blue
             1.00 0.55 0.10     % orange
             0.85 0.30 0.90     % magenta
             1.00 1.00 0.30     % yellow
             0.30 1.00 0.80     % teal
             1.00 0.40 0.40     % salmon
             0.60 0.60 1.00 ];  % periwinkle
    if n <= size(base,1), cmap = base(1:n,:);
    else, cmap = base(mod(0:n-1, size(base,1))+1, :); end
end

function [bx, by] = bed_endpoints(calib, W, H)
    A = calib.lineA; B = calib.lineB; C = calib.lineC;
    if abs(B) > eps                       % non-vertical: span full width
        bx = [1 W];
        by = -(A*bx + C) / B;
    else                                  % vertical line: x = -C/A
        xv = -C / A;
        bx = [xv xv]; by = [1 H];
    end
end

function img = draw_dashed_line(img, bx, by, col, thick)
    H = size(img,1); W = size(img,2);
    n = max(abs(diff(bx)), abs(diff(by))); n = max(round(n),1);
    xs = linspace(bx(1), bx(2), n+1);
    ys = linspace(by(1), by(2), n+1);
    dash = 10*thick; gap = 7*thick;        % dash pattern in px
    for i = 1:numel(xs)
        if mod(i-1, dash+gap) < dash
            img = stamp(img, xs(i), ys(i), thick, col, W, H);
        end
    end
end

function img = draw_circle(img, cx, cy, rad, col, thick)
    H = size(img,1); W = size(img,2);
    th = linspace(0, 2*pi, max(24, round(2*pi*rad)));
    for t = th
        img = stamp(img, cx + rad*cos(t), cy + rad*sin(t), thick, col, W, H);
    end
end

function img = draw_disc(img, cx, cy, rad, col)
    H = size(img,1); W = size(img,2);
    x0 = max(1,floor(cx-rad)); x1 = min(W,ceil(cx+rad));
    y0 = max(1,floor(cy-rad)); y1 = min(H,ceil(cy+rad));
    for x = x0:x1
        for y = y0:y1
            if (x-cx)^2 + (y-cy)^2 <= rad^2
                img(y,x,1)=col(1); img(y,x,2)=col(2); img(y,x,3)=col(3);
            end
        end
    end
end

function img = draw_border(img, col, thick)
    img(1:thick,:,1)=col(1);        img(1:thick,:,2)=col(2);        img(1:thick,:,3)=col(3);
    img(end-thick+1:end,:,1)=col(1);img(end-thick+1:end,:,2)=col(2);img(end-thick+1:end,:,3)=col(3);
    img(:,1:thick,1)=col(1);        img(:,1:thick,2)=col(2);        img(:,1:thick,3)=col(3);
    img(:,end-thick+1:end,1)=col(1);img(:,end-thick+1:end,2)=col(2);img(:,end-thick+1:end,3)=col(3);
end

function img = stamp(img, x, y, thick, col, W, H)
    xr = max(1,round(x)-thick+1):min(W,round(x)+thick-1);
    yr = max(1,round(y)-thick+1):min(H,round(y)+thick-1);
    if isempty(xr) || isempty(yr), return; end
    img(yr,xr,1)=col(1); img(yr,xr,2)=col(2); img(yr,xr,3)=col(3);
end

function img = local_read_window_frame(vr, wStart, k, filterType)
%LOCAL_READ_WINDOW_FRAME  Window index k -> the matching frame of the video.
%   Window index k is absolute video frame wStart + k - 1. Seeks per call,
%   which is fine here: QA renders at most a few hundred frames, unlike the
%   detection pass which streams sequentially.
    absFrame      = wStart + k - 1;
    vr.CurrentTime = (absFrame - 1) / vr.FrameRate;
    if ~hasFrame(vr)
        error('make_annotated_frames:frameBeyondEnd', ...
              'Video frame %d (window index %d) is past the end.', absFrame, k);
    end
    img = readFrame(vr);
    if nargin >= 4 && ~isempty(filterType)
        img = apply_filter(img, filterType);
    end
end
