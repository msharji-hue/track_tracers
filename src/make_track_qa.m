function make_track_qa(meta, tracks, nDetected, outRoot)
% MAKE_TRACK_QA  Compact static QA image per trial (saved off-screen).
%
%   (1) tracked marker x vs tracking frame   (2) tracked marker y vs tracking frame
%   (3) detections per ORIGINAL video frame  (4) trial metadata summary
%
%   Tracking-only: no kinematics, no bending panels.
%
%   Indexing: tracks.trackedX columns are TRACKING-frame index (column 1 ==
%   firstValidFrame). nDetected is indexed by EXPORTED-WINDOW frame, i.e. the
%   sorted PNG list in framesDir. Both are window-relative, so they plot
%   against each other directly; neither is an absolute video-frame index. Add
%   meta.windowStart - 1 to reach the video frame.
%
%   NOTE: uses a fixed RGB palette rather than lines(), which is a colormap
%   function unavailable in some headless/-batch MATLAB sessions (this caused
%   "Unrecognized function or variable 'lines'").

    qaDir = fullfile(outRoot, 'qa');
    if ~exist(qaDir, 'dir'), mkdir(qaDir); end

    f = figure('Visible','off', 'Name','Track QA', 'Position',[100 100 1250 760]);
    nM   = size(tracks.trackedX, 1);
    cmap = marker_palette(nM);

    % ── (1) marker x vs tracking frame ───────────────────────────────────
    subplot(2,2,1); hold on
    for m = 1:nM
        plot(tracks.trackedX(m,:), '.-', 'Color', cmap(m,:), 'MarkerSize',4);
    end
    xline(1, 'g', 'trackFrame 1', 'LabelVerticalAlignment','bottom');
    grid on; xlabel('tracking frame (1 = firstValidFrame)'); ylabel('x (px)');
    title(sprintf('%s   —   firstValidFrame = %d', ...
        strrep(meta.trialTag,'_','\_'), meta.firstValidFrame));

    % ── (2) marker y vs tracking frame ───────────────────────────────────
    subplot(2,2,2); hold on
    for m = 1:nM
        plot(tracks.trackedY(m,:), '.-', 'Color', cmap(m,:), 'MarkerSize',4);
    end
    grid on; xlabel('tracking frame'); ylabel('y (px)');
    title('tracked marker y (per-ID colour: watch for ID swaps)');

    % ── (3) detections per original frame ────────────────────────────────
    subplot(2,2,3);
    bar(nDetected, 'FaceColor',[.30 .60 .90], 'EdgeColor','none');
    yline(8, 'k--', '8 markers');
    if isfinite(meta.firstValidFrame)
        xline(meta.firstValidFrame, 'g', 'firstValid', 'LabelVerticalAlignment','bottom');
    end
    grid on; xlabel('original video frame'); ylabel('n detected');
    title('detections per frame');

    % ── (4) metadata summary ─────────────────────────────────────────────
    subplot(2,2,4); axis off
    txt = {
        sprintf('material        : %s', meta.material)
        sprintf('condition       : %s', meta.container)
        sprintf('batch           : %s', meta.batchName)
        sprintf('drop height     : %g mm', meta.dropHeight_mm)
        sprintf('trial           : T%02d', meta.trialNum)
        sprintf('fps             : %.2f', meta.fps_true)
        sprintf('nFrames (video) : %s', num2str(getf(meta,'nFrames')))
        sprintf('firstValidFrame : %d', meta.firstValidFrame)
        sprintf('nTracked frames : %s', num2str(getf(meta,'nTracked')))
        sprintf('markers         : %d', nM)
        sprintf('rho_particle    : %.2f g/cm^3', getf(meta,'rho_particle_g_cm3'))
        sprintf('rho_bulk        : %.3f g/cm^3', getf(meta,'rho_bulk_g_cm3'))
        sprintf('phi             : %.3f', getf(meta,'phi')) };
    text(0.02, 0.98, txt, 'VerticalAlignment','top', ...
        'FontName','FixedWidth', 'FontSize',10, 'Interpreter','none');
    title('trial metadata');

    outPng = fullfile(qaDir, [meta.trialTag '_qa.png']);
    exportgraphics(f, outPng, 'Resolution', 150);
    close(f);
    fprintf('  Saved QA figure  : %s\n', outPng);
end

function cmap = marker_palette(n)
% Fixed palette — no dependence on lines()/colormap functions.
    base = [ 0.20 0.90 0.20
             0.20 0.70 1.00
             1.00 0.55 0.10
             0.85 0.30 0.90
             0.90 0.90 0.20
             0.30 0.85 0.75
             1.00 0.40 0.40
             0.60 0.60 1.00 ];
    if n <= size(base,1), cmap = base(1:n,:);
    else, cmap = base(mod(0:n-1, size(base,1))+1, :); end
end

function v = getf(s, f)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end
