function make_annotated_video(framesDir, det, kin, params, saveVideo)
% MAKE_ANNOTATED_VIDEO  Annotated video with side-by-side kinematic plots.
%
%   Layout: video (left, rotated 90° CCW) | z(t) and v(t) plots (right)
%   White background, publishable style.
%
%   Inputs:
%       framesDir  - path to folder containing frame_XXXXX.png files
%       det        - detections struct
%       kin        - kinematics struct
%       params     - parameter struct (built by scripts/make_video.m)
%       saveVideo  - true to save .mp4, false to preview only

    % ── Setup ─────────────────────────────────────────────────────────────
    files      = sort({dir(fullfile(framesDir,'*.png')).name});
    fps_true   = params.fps_true;
    outputFps  = params.outputFps;
    slowFactor = params.slowFactor;
    startFrame = params.startFrame;
    stopFrame  = min(params.stopFrame, numel(files));
    impactIdx  = kin.impact_index;
    mmPerPx    = params.mmPerPx;
    trailLen   = params.trailLen;
    trailSize  = params.trailSize;
    nFrames    = stopFrame - startFrame + 1;
    nKin       = numel(kin.t_s);

    effectiveSlowdown = fps_true / outputFps * slowFactor;
    fprintf('Effective slowdown: %.0fx  (fps_true=%.0f, outputFps=%d, slowFactor=%d)\n',...
        effectiveSlowdown, fps_true, outputFps, slowFactor);

    % Original frame size
    tmp = imread(fullfile(framesDir, files{startFrame}));
    fH  = size(tmp,1);
    fW  = size(tmp,2);

    % After CCW rotation: fW rows × fH cols
    rfH = fW;
    rfW = fH;

    % ── Display geometry ──────────────────────────────────────────────────
    vid_disp_H = 520;
    sy         = vid_disp_H / rfH;
    content_W  = round(rfW * sy);
    vid_disp_W = max(180, content_W + 50);
    pad_left   = floor((vid_disp_W - content_W) / 2);
    pan_W      = 720;
    figW       = vid_disp_W + pan_W;
    figH       = vid_disp_H;

    % Scale bar (10 mm)
    sbPx = round(10 / mmPerPx * sy);

    % ── Coordinate transforms (CCW, YDir='reverse') ───────────────────────
    orig2x = @(xo, yo)  yo .* sy + pad_left;
    orig2y = @(xo,  ~) (rfH - xo) .* sy;

    % ── Figure ────────────────────────────────────────────────────────────
    hFig = figure('Name','Trial Video','Units','pixels',...
                  'Position',[50 80 figW figH],'Color','w',...
                  'ToolBar','none','MenuBar','none');

    % ── Video axes ────────────────────────────────────────────────────────
    ax_vid = axes(hFig,'Units','pixels','Position',[0 0 vid_disp_W vid_disp_H]);
    set(ax_vid,'Color','w','XColor','none','YColor','none','Visible','off');
    hold(ax_vid,'on');

    % ── Plot axes ─────────────────────────────────────────────────────────
    ml=68; mr=18; mb=42; mt=14; gap=22;
    ph = floor((figH - mt - mb - gap) / 2);
    pw = pan_W - ml - mr;
    ax_z = axes(hFig,'Units','pixels','Position',[vid_disp_W+ml, mb+ph+gap, pw, ph]);
    ax_v = axes(hFig,'Units','pixels','Position',[vid_disp_W+ml, mb,         pw, ph]);

    pStyle = {'Color',      'white', ...
              'XColor',     [0.10 0.10 0.10], ...
              'YColor',     [0.10 0.10 0.10], ...
              'FontSize',   12, ...
              'LineWidth',   1.2, ...
              'Box',        'on', ...
              'GridAlpha',   0.13, ...
              'GridColor',  [0.30 0.30 0.30], ...
              'TickDir',    'in', ...
              'XMinorTick', 'on', ...
              'YMinorTick', 'on'};
    set(ax_z, pStyle{:}); hold(ax_z,'on'); grid(ax_z,'on');
    set(ax_v, pStyle{:}); hold(ax_v,'on'); grid(ax_v,'on');

    % ── Kinematics (clip at kin.stopFrame — post-stop data is NaN) ────────
    preFrames   = round(fps_true * 0.008);
    plotStart   = max(1, impactIdx - preFrames);
    plotEnd     = min(kin.stopFrame, nKin);
    tKin_ms     = ((1:nKin) - impactIdx) ./ fps_true .* 1000;
    % kd_kinematics saves the smoothed velocity as kin.v and the depth as
    % kin.z (an alias for kin.depthRod_cm). v_smooth/z_smooth are internal
    % variable names inside kd_kinematics and were never struct fields, so
    % reading them here returned nothing.
    vAll        = kin.v;
    zAll        = kin.z;
    t_sl        = tKin_ms(plotStart:plotEnd);
    z_sl        = zAll(plotStart:plotEnd);
    v_sl        = vAll(plotStart:plotEnd);
    t_stop_ms   = tKin_ms(plotEnd);

    % Static black traces
    plot(ax_z, t_sl, z_sl, 'k-', 'LineWidth',1.5);
    plot(ax_v, t_sl, v_sl, 'k-', 'LineWidth',1.5);

    % Impact: grey dashed (label on z only)
    xline(ax_z, 0, '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1.0, ...
        'Label','impact', 'LabelVerticalAlignment','top', ...
        'LabelHorizontalAlignment','right', 'FontSize',9);
    xline(ax_v, 0, '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1.0);

    % Stop: red dashed (label on z only)
    xline(ax_z, t_stop_ms, '--', 'Color',[0.80 0.10 0.10], 'LineWidth',1.0, ...
        'Label','stop', 'LabelVerticalAlignment','top', ...
        'LabelHorizontalAlignment','left', 'FontSize',9);
    xline(ax_v, t_stop_ms, '--', 'Color',[0.80 0.10 0.10], 'LineWidth',1.0);

    % Zero velocity reference
    yline(ax_v, 0, '--', 'Color',[0.65 0.65 0.65], 'LineWidth',0.8);

    % Axis limits — clip a few ms after stop
    clipMs = 5;
    xlim(ax_z, [t_sl(1), t_stop_ms + clipMs]);
    xlim(ax_v, [t_sl(1), t_stop_ms + clipMs]);
    zpad = max(0.06*range(z_sl), 0.05);
    ylim(ax_z, [min(z_sl)-zpad,  max(z_sl)+zpad]);
    vpad = max(0.10*range(v_sl), 1);
    ylim(ax_v, [min(v_sl)-vpad,  max(v_sl)+vpad]);

    % Labels — plain TeX
    set(ax_z,'XTickLabel',{});
    ylabel(ax_z, 'z  (cm)',         'FontSize',13, 'Interpreter','tex');
    xlabel(ax_v, 't  (ms)',         'FontSize',13, 'Interpreter','tex');
    ylabel(ax_v, 'v  (cm s^{-1})', 'FontSize',13, 'Interpreter','tex');

    % Moving dots
    hDot_z = plot(ax_z, nan, nan, 'o', 'MarkerFaceColor',[0.15 0.40 0.80],...
        'MarkerEdgeColor','none', 'MarkerSize',7);
    hDot_v = plot(ax_v, nan, nan, 'o', 'MarkerFaceColor',[0.15 0.40 0.80],...
        'MarkerEdgeColor','none', 'MarkerSize',7);

    % Readout boxes (top-left to avoid crowding stop line)
    hTxt_z = text(ax_z, 0.03, 0.92, 'z = 0.00 cm', 'Units','normalized',...
        'HorizontalAlignment','left', 'FontSize',11, 'FontWeight','bold',...
        'Color',[0.10 0.10 0.10], 'BackgroundColor',[1 1 1 0.85],...
        'EdgeColor',[0.50 0.50 0.50], 'Margin',3, 'Interpreter','none');
    hTxt_v = text(ax_v, 0.03, 0.92, 'v = 0.0 cm/s', 'Units','normalized',...
        'HorizontalAlignment','left', 'FontSize',11, 'FontWeight','bold',...
        'Color',[0.10 0.10 0.10], 'BackgroundColor',[1 1 1 0.85],...
        'EdgeColor',[0.50 0.50 0.50], 'Margin',3, 'Interpreter','none');

    % ── Static video overlays ─────────────────────────────────────────────
    rgb0     = rot90(imread(fullfile(framesDir,files{startFrame})), 1);
    rgb0_s   = imresize(rgb0, [vid_disp_H NaN]);
    rgb0_pad = pad_lr(rgb0_s, vid_disp_W, pad_left);
    hImg     = image(ax_vid, rgb0_pad);
    set(ax_vid,'XLim',[0 vid_disp_W],'YLim',[0 vid_disp_H],'YDir','reverse');

    % Scale bar
    sbX2 = vid_disp_W - pad_left - 8;
    sbX1 = sbX2 - sbPx;
    sbY  = vid_disp_H - 16;
    plot(ax_vid,[sbX1 sbX2],[sbY sbY],   '-','Color','w','LineWidth',2.2);
    plot(ax_vid,[sbX1 sbX1],[sbY-4 sbY+4],'-','Color','w','LineWidth',1.5);
    plot(ax_vid,[sbX2 sbX2],[sbY-4 sbY+4],'-','Color','w','LineWidth',1.5);
    text(ax_vid,(sbX1+sbX2)/2, sbY-10, '10 mm', 'Color','w', 'FontSize',8,...
        'FontWeight','bold', 'HorizontalAlignment','center');

    % Bed line
    plot(ax_vid,...
        [orig2x(params.bedPoint1(1),params.bedPoint1(2)), ...
         orig2x(params.bedPoint2(1),params.bedPoint2(2))], ...
        [orig2y(params.bedPoint1(1),params.bedPoint1(2)), ...
         orig2y(params.bedPoint2(1),params.bedPoint2(2))], ...
        '--','Color',params.bedColor,'LineWidth',2.0);

    % Time stamp
    hTxt_t = text(ax_vid, pad_left+4, 15, 't = 0.0 ms', ...
        'Color','white', 'FontSize',10, 'FontWeight','bold',...
        'BackgroundColor',[0 0 0 0.55]);

    % ── VideoWriter ───────────────────────────────────────────────────────
    if saveVideo
        vw           = VideoWriter(params.outputName,'MPEG-4');
        vw.FrameRate = outputFps;
        vw.Quality   = 95;
        open(vw);
        fprintf('Saving: %d frames × %d repeats at %d fps = %.1f s video\n',...
            nFrames, slowFactor, outputFps, nFrames*slowFactor/outputFps);
    end

    theta     = linspace(0,2*pi,60);
    pauseTime = 1/outputFps;
    hDynamic  = gobjects(0);

    % ── Render loop ───────────────────────────────────────────────────────
    for i = startFrame:stopFrame
        if ~ishandle(hFig), break; end

        rgb     = rot90(imread(fullfile(framesDir,files{i})), 1);
        rgb_s   = imresize(rgb, [vid_disp_H NaN]);
        rgb_pad = pad_lr(rgb_s, vid_disp_W, pad_left);
        centers = det.detect.centersCell{i};
        radii   = det.detect.radiiCell{i};
        t_ms    = (i - impactIdx) / fps_true * 1000;

        % Fast image update
        set(hImg,'CData',rgb_pad);

        % Remove previous dynamic overlays
        delete(hDynamic(ishandle(hDynamic)));
        hDynamic = gobjects(0);

        % Impact flash
        if i == impactIdx
            h = plot(ax_vid,[1 vid_disp_W vid_disp_W 1 1],...
                            [1 1 vid_disp_H vid_disp_H 1],...
                '-','Color',[1 0.30 0],'LineWidth',7);
            hDynamic(end+1) = h;
        end

        % Circles + fading trails
        if ~isempty(centers)
            for ci = 1:size(centers,1)
                cxo = centers(ci,1);  cyo = centers(ci,2);
                cxd = orig2x(cxo,cyo);
                cyd = orig2y(cxo,cyo);
                rd  = radii(ci) * sy;

                h = plot(ax_vid, cxd+rd.*cos(theta), cyd+rd.*sin(theta), '-',...
                    'Color',params.markerColor,'LineWidth',params.circLW);
                hDynamic(end+1) = h;

                for ti = max(startFrame,i-trailLen):i-1
                    cp = det.detect.centersCell{ti};
                    if isempty(cp)||ci>size(cp,1), continue; end
                    a = params.trailAlpha*(ti-max(startFrame,i-trailLen)+1)/trailLen;
                    h = plot(ax_vid,...
                        orig2x(cp(ci,1),cp(ci,2)),...
                        orig2y(cp(ci,1),cp(ci,2)),'o',...
                        'MarkerSize',      trailSize,...
                        'MarkerFaceColor', params.markerColor*a,...
                        'MarkerEdgeColor', 'none');
                    hDynamic(end+1) = h;
                end
            end
        end

        % Update plots — clamp to valid kinematics range
        i_kin = min(i, plotEnd);
        set(hTxt_t,'String', sprintf('t = %.1f ms', t_ms));
        set(hDot_z,'XData',  tKin_ms(i_kin), 'YData', zAll(i_kin));
        set(hDot_v,'XData',  tKin_ms(i_kin), 'YData', vAll(i_kin));
        set(hTxt_z,'String', sprintf('z = %.2f cm',  zAll(i_kin)));
        set(hTxt_v,'String', sprintf('v = %.1f cm/s', max(0, vAll(i_kin))));

        drawnow;

        if saveVideo
            F = getframe(hFig);
            for rep = 1:slowFactor
                writeVideo(vw, F);
            end
            if mod(i-startFrame,50)==0||i==stopFrame
                fprintf('  %d/%d  (%.0f%%)\n',...
                    i-startFrame+1,nFrames,100*(i-startFrame+1)/nFrames);
            end
        else
            pause(pauseTime);
        end
    end

    if saveVideo
        close(vw);
        fprintf('Done: %s\n', params.outputName);
        fprintf('Video duration: %.1f s  (%.0fx slowdown)\n',...
            nFrames*slowFactor/outputFps, effectiveSlowdown);
    end
    if ishandle(hFig)&&~saveVideo, close(hFig); end
end

% ── Helper: pad image left/right with white ────────────────────────────────
function out = pad_lr(img, panel_W, pad_left)
    [ih, iw, ch] = size(img);
    pad_r = max(0, panel_W - pad_left - iw);
    out   = [ones(ih,pad_left,ch,'uint8')*255, img, ...
             ones(ih,pad_r,   ch,'uint8')*255];
end