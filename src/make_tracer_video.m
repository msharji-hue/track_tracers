function make_tracer_video(framesDir, detectOut, params, doSave)
% MAKE_TRACER_VIDEO  Qualitative display clip: rotated video with tracer
% circles + dashed bed line, and a bold TWO-LINE info panel above the video.
% No plots, no kinematics.
%
%   Everything (panel + video + overlays) lives in ONE axes that fills the
%   figure, captured with a single getframe -> the panel can never be cropped.
%
%   Because the rotated frame can be very narrow, the canvas is widened to
%   params.minWidth and the video is centered, so the caption has room to be
%   readable. Caption is two lines, auto-sized to fit the panel.
%
%   Orientation matches make_annotated_video.m: each frame is rot90(img,1), with
%       x_rot = y_orig ,  y_rot = (rotatedHeight - x_orig)
%
%   params fields:
%       outputFps, slowFactor, fps_true
%       markerColor, bedColor, circLW
%       bedPoint1, bedPoint2   [x y] in ORIGINAL (unrotated) pixel coords
%       trailLen, trailAlpha, trailSize
%       rotateCCW (default true)
%       panelColor (also the side-bar colour)
%       panelFrac  panel height as fraction of frame height (default 0.10)
%       minWidth   minimum canvas width in px (default 220)
%       outputName full path to .mp4 (doSave = true)

    if nargin < 4, doSave = false; end
    if ~isfield(params,'markerColor'), params.markerColor = [0.35 0.95 0.25]; end
    if ~isfield(params,'bedColor'),    params.bedColor    = [1.00 0.90 0.20]; end
    if ~isfield(params,'circLW'),      params.circLW      = 2;                end
    if ~isfield(params,'trailLen'),    params.trailLen    = 8;                end
    if ~isfield(params,'trailAlpha'),  params.trailAlpha  = 0.5;              end
    if ~isfield(params,'trailSize'),   params.trailSize   = 6;                end
    if ~isfield(params,'rotateCCW'),   params.rotateCCW   = true;             end
    if ~isfield(params,'panelColor'),  params.panelColor  = [0.10 0.10 0.10]; end
    if ~isfield(params,'panelFrac'),   params.panelFrac   = 0.10;             end
    if ~isfield(params,'minWidth'),    params.minWidth    = 220;              end
    if ~isfield(params,'fps_true'),    params.fps_true    = params.outputFps; end

    files   = dir(fullfile(framesDir, '*.png'));
    files   = sort({files.name});
    nFrames = min(numel(files), numel(detectOut.centersCell));
    if nFrames == 0, error('No frames / detections to render.'); end

    % ── geometry ─────────────────────────────────────────────────────────
    img0    = readRot(fullfile(framesDir, files{1}), params.rotateCCW);
    rfH     = size(img0,1);   rfW = size(img0,2);
    panelH  = max(64, round(params.panelFrac * rfH));
    canvasW = max(rfW, params.minWidth);          % widen so caption fits
    xPad    = floor((canvasW - rfW) / 2);         % centre the video
    totalH  = rfH + panelH;

    ss    = get(0,'ScreenSize');
    scale = min([1, 0.9*ss(3)/canvasW, 0.85*ss(4)/totalH]);
    figW  = max(1, round(canvasW * scale));
    figH  = max(1, round(totalH  * scale));

    fig = figure('Visible','on', 'Color',params.panelColor, 'MenuBar','none', ...
                 'ToolBar','none', 'Name','Tracer video', 'Resize','off', ...
                 'Position',[60 60 figW figH]);
    ax = axes('Parent',fig, 'Units','normalized', 'Position',[0 0 1 1]);

    % caption font (screen px) sized to fit two lines inside the panel
    sample = {sprintf('t = %.2f ms', 88.88), sprintf('frame %d/%d', nFrames, nFrames)};
    Nmax   = max(cellfun(@numel, sample));
    fByW   = (figW * 0.90) / (0.62 * Nmax);       % width limit (longest line)
    fByH   = (panelH * scale * 0.82) / (2 * 1.30); % height limit (2 lines)
    fontPx = max(9, floor(min(fByW, fByH)));

    if doSave
        vw = VideoWriter(params.outputName, 'MPEG-4');
        vw.FrameRate = params.outputFps; vw.Quality = 95;
        open(vw);
    end

    % bed line: rotate endpoints, offset into canvas, extend across frame
    [bx, by] = orig2rot(params.bedPoint1(1), params.bedPoint1(2), rfH, params.rotateCCW);
    [cx, cy] = orig2rot(params.bedPoint2(1), params.bedPoint2(2), rfH, params.rotateCCW);
    bx = bx + xPad; cx = cx + xPad; by = by + panelH; cy = cy + panelH;
    d = [cx-bx, cy-by]; if norm(d) > 0, d = d/norm(d); end
    bedX = [bx - 5000*d(1), cx + 5000*d(1)];
    bedY = [by - 5000*d(2), cy + 5000*d(2)];

    targetSize = [];
    for i = 1:nFrames
        img = readRot(fullfile(framesDir, files{i}), params.rotateCCW);

        cla(ax);
        % video in the lower band, centred; top + sides show panelColor
        image(ax, [xPad+1 xPad+rfW], [panelH+1 totalH], img);
        set(ax, 'YDir','reverse', 'XLim',[0.5 canvasW+0.5], 'YLim',[0.5 totalH+0.5], ...
                'Color',params.panelColor, 'XColor','none', 'YColor','none', ...
                'XTick',[], 'YTick',[], 'Position',[0 0 1 1]);
        hold(ax,'on');

        plot(ax, bedX, bedY, '--', 'Color',params.bedColor, 'LineWidth',2);

        if params.trailLen > 0
            tx = []; ty = [];
            for j = max(1, i-params.trailLen):(i-1)
                c = detectOut.centersCell{j};
                if ~isempty(c)
                    [rx, ry] = orig2rot(c(:,1), c(:,2), rfH, params.rotateCCW);
                    tx = [tx; rx+xPad]; ty = [ty; ry+panelH]; %#ok<AGROW>
                end
            end
            if ~isempty(tx)
                scatter(ax, tx, ty, params.trailSize, 'filled', ...
                    'MarkerFaceColor',params.markerColor, ...
                    'MarkerFaceAlpha',params.trailAlpha, 'MarkerEdgeColor','none');
            end
        end

        c = detectOut.centersCell{i};
        r = detectOut.radiiCell{i};
        if ~isempty(c)
            [rx, ry] = orig2rot(c(:,1), c(:,2), rfH, params.rotateCCW);
            viscircles(ax, [rx+xPad ry+panelH], r, 'Color',params.markerColor, ...
                'LineWidth',params.circLW, 'EnhanceVisibility',false);
            plot(ax, rx+xPad, ry+panelH, '.', 'Color',params.markerColor, 'MarkerSize',8);
        end

        % two-line caption, centred in the panel band
        t_ms = (i-1) / params.fps_true * 1000;
        cap  = {sprintf('t = %.2f ms', t_ms); sprintf('frame %d/%d', i, nFrames)};
        text(ax, canvasW/2, panelH/2, cap, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'Color','w', 'FontWeight','bold', 'FontUnits','pixels', ...
            'FontSize',fontPx, 'Clipping','off', 'Interpreter','none');
        hold(ax,'off'); drawnow;

        if doSave
            F = getframe(ax); frame = F.cdata;
            if isempty(targetSize)
                targetSize = [size(frame,1) size(frame,2)];
            elseif ~isequal([size(frame,1) size(frame,2)], targetSize)
                frame = imresize(frame, targetSize);
            end
            for k = 1:params.slowFactor
                writeVideo(vw, frame);
            end
        else
            pause(1 / params.outputFps);
        end
    end

    if doSave
        close(vw);
        fprintf('Saved video: %s\n', params.outputName);
    end
    if ishandle(fig), close(fig); end
end

% ── helpers ───────────────────────────────────────────────────────────────
function img = readRot(fpath, rotateCCW)
    img = imread(fpath);
    if rotateCCW, img = rot90(img, 1); end
end

function [xr, yr] = orig2rot(xo, yo, rfH, rotateCCW)
    if rotateCCW
        xr = yo;
        yr = rfH - xo;
    else
        xr = xo;
        yr = yo;
    end
end