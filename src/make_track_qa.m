function make_track_qa(meta, tracks, nDetected, outRoot)
% MAKE_TRACK_QA  One compact QA image per trial (saved off-screen).
%
%   2 x 3 layout:
%     (1) marker x vs tracking frame        (4) detections per original frame
%     (2) bending amplitude (mm) vs frame   (5) rod shape: reference vs peak
%     (3) bending/tilt angle (deg) vs frame (6) peak-value summary
%
%   Indexing: tracks.trackedX columns are TRACKING-frame index (col 1 ==
%   firstValidFrame). nDetected is ORIGINAL video-frame index.

    qaDir = fullfile(outRoot, 'qa');
    if ~exist(qaDir, 'dir'), mkdir(qaDir); end

    f = figure('Visible','off', 'Name','Track QA', 'Position',[80 80 1500 780]);
    nM   = size(tracks.trackedX, 1);
    cmap = lines(max(nM, 1));
    hasBend = isfield(tracks,'bending') && ~isfield(tracks.bending,'error') ...
              && isfield(tracks.bending,'rms_mm');
    if hasBend, b = tracks.bending; end

    % ── (1) marker x vs tracking frame ───────────────────────────────────
    subplot(2,3,1); hold on
    for m = 1:nM
        plot(tracks.trackedX(m,:), '.-', 'Color', cmap(m,:));
    end
    xline(1, 'g', 'firstValid', 'LabelVerticalAlignment','bottom');
    if isfinite(tracks.impact_index)
        xline(tracks.impact_index, '--', 'Color',[.5 .5 .5], 'Label','impact', ...
              'LabelVerticalAlignment','bottom');
    end
    if isfinite(tracks.stopFrame)
        xline(tracks.stopFrame, 'r', 'Label','stop', 'LabelVerticalAlignment','bottom');
    end
    grid on; xlabel('tracking frame'); ylabel('x (px)');
    title(sprintf('%s   fvf=%d', strrep(meta.trialTag,'_','\_'), meta.firstValidFrame));

    % ── (2) bending amplitude (mm) vs frame ──────────────────────────────
    subplot(2,3,2); hold on
    if hasBend
        plot(b.rms_mm, '-', 'Color',[0.00 0.45 0.74], 'DisplayName','RMS');
        plot(b.max_mm, '-', 'Color',[0.85 0.33 0.10], 'DisplayName','max');
        if isfinite(b.baseline_rms_mm)
            yline(b.baseline_rms_mm, ':', 'baseline', 'Color',[.3 .3 .3]);
        end
        markLines(tracks, b);
        legend('Location','best'); grid on
        xlabel('tracking frame'); ylabel('deflection (mm)');
        title(sprintf('bending mm — peak RMS %.3f (flag %d)', b.peak_rms_mm, b.bendFlag));
    else, noData(); end

    % ── (3) angles (deg) vs frame ────────────────────────────────────────
    subplot(2,3,3); hold on
    if hasBend
        plot(b.tilt_deg,       '-', 'Color',[0.49 0.18 0.56], 'DisplayName','tilt');
        plot(b.bend_angle_deg, '-', 'Color',[0.00 0.50 0.20], 'DisplayName','bend angle');
        plot(b.seg_angle_deg,  ':', 'Color',[0.6 0.6 0.2],   'DisplayName','seg max');
        markLines(tracks, b);
        legend('Location','best'); grid on
        xlabel('tracking frame'); ylabel('angle (deg)');
        title(sprintf('angles — tilt %.2f, bend %.2f deg (tiltflag %d)', ...
              b.tilt_peak_deg, b.bend_angle_peak_deg, b.tiltFlag));
    else, noData(); end

    % ── (4) detections per original frame ────────────────────────────────
    subplot(2,3,4);
    bar(nDetected, 'FaceColor',[.30 .60 .90], 'EdgeColor','none');
    yline(8, 'k--', '8 markers');
    if isfinite(meta.firstValidFrame)
        xline(meta.firstValidFrame, 'g', 'firstValid', 'LabelVerticalAlignment','bottom');
    end
    grid on; xlabel('original video frame'); ylabel('n detected'); title('detections per frame');

    % ── (5) rod shape: reference vs peak ─────────────────────────────────
    subplot(2,3,5); hold on
    if hasBend
        rf = b.refFrame; pf = b.peakFrame;
        if isfinite(rf)
            plot(tracks.trackedX(:,rf), tracks.trackedY(:,rf), 'o-', ...
                 'Color',[.45 .45 .45], 'DisplayName','reference');
        end
        if isfinite(pf)
            plot(tracks.trackedX(:,pf), tracks.trackedY(:,pf), 'o-', ...
                 'Color',[0.85 0.10 0.10], 'DisplayName','peak bend');
        end
        set(gca,'YDir','reverse'); axis equal; grid on
        xlabel('x (px)'); ylabel('y (px)'); legend('Location','best');
        title('rod shape: reference vs peak');
    else, noData(); end

    % ── (6) peak-value summary ───────────────────────────────────────────
    subplot(2,3,6); axis off
    if hasBend
        lines = {
            sprintf('peak RMS deflection : %.3f mm', b.peak_rms_mm)
            sprintf('peak max deflection : %.3f mm', b.peak_max_mm)
            sprintf('bend angle (peak)   : %.2f deg', b.bend_angle_peak_deg)
            sprintf('seg-angle max       : %.2f deg', b.seg_angle_peak_deg)
            sprintf('rod tilt (peak)     : %.2f deg', b.tilt_peak_deg)
            sprintf('curvature (peak)    : %.4g 1/mm', b.curv_peak_1pmm)
            sprintf('bend @ stop         : %.3f mm', b.bend_at_stop_mm)
            sprintf('baseline RMS        : %.3f mm', b.baseline_rms_mm)
            sprintf('t(peak) - impact    : %.2f ms', b.t_peak_ms)
            sprintf('bendFLAG=%d   tiltFLAG=%d', b.bendFlag, b.tiltFlag) };
        text(0.02, 0.98, lines, 'VerticalAlignment','top', ...
             'FontName','FixedWidth', 'FontSize',10, 'Interpreter','none');
    end
    title('bending summary');

    outPng = fullfile(qaDir, [meta.trialTag '_qa.png']);
    exportgraphics(f, outPng, 'Resolution', 150);
    close(f);
    fprintf('  Saved QA figure  : %s\n', outPng);
end

function markLines(tracks, b)
    if isfinite(tracks.impact_index), xline(tracks.impact_index, '--', 'Color',[.5 .5 .5]); end
    if isfinite(tracks.stopFrame),    xline(tracks.stopFrame, 'r'); end
    if isfinite(b.peakFrame),         xline(b.peakFrame, 'm', 'peak'); end
end

function noData()
    text(0.5,0.5,'no bending data','HorizontalAlignment','center'); axis off
end
