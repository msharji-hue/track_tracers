function endFrame = pick_end_frame(v, startFrame)
% PICK_END_FRAME  Visually browse frames and confirm your end frame.
%   Same controls as pick_start_frame but cannot go below startFrame.

    totalFrames = floor(v.Duration * v.FrameRate);
    endFrame    = totalFrames;  % default to last frame
    confirmed   = false;

    % ── Build figure ─────────────────────────────────────────
    fig = figure('Name', 'Pick End Frame', ...
                 'NumberTitle', 'off', ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'none', ...
                 'Position', [100 100 800 600]);

    ax = axes('Parent', fig, 'Position', [0.05 0.2 0.9 0.75]);

    % ── Buttons ──────────────────────────────────────────────
    uicontrol('Style', 'pushbutton', 'String', '<<',  ...
        'Position', [30  20 60 40],  'Callback', @(~,~) step(-10));
    uicontrol('Style', 'pushbutton', 'String', '<',   ...
        'Position', [100 20 60 40],  'Callback', @(~,~) step(-1));
    uicontrol('Style', 'pushbutton', 'String', '>',   ...
        'Position', [170 20 60 40],  'Callback', @(~,~) step(1));
    uicontrol('Style', 'pushbutton', 'String', '>>',  ...
        'Position', [240 20 60 40],  'Callback', @(~,~) step(10));
    uicontrol('Style', 'pushbutton', 'String', 'Jump to Frame', ...
        'Position', [320 20 120 40], 'Callback', @(~,~) jump());
    uicontrol('Style', 'pushbutton', 'String', '✓ Confirm End Frame', ...
        'Position', [460 20 200 40], ...
        'BackgroundColor', [0.8 0.2 0.2], 'ForegroundColor', 'white', ...
        'FontWeight', 'bold', ...
        'Callback', @(~,~) confirm());

    show_frame();
    waitfor(fig);

    if ~confirmed
        error('No end frame selected — figure was closed without confirming.');
    end

    % ── Nested helpers ────────────────────────────────────────
    function show_frame()
        v.CurrentTime = (endFrame - 1) / v.FrameRate;
        imshow(readFrame(v), 'Parent', ax);
        title(ax, sprintf('Frame  %d / %d    (start: %d)', endFrame, totalFrames, startFrame), ...
              'FontSize', 14);
        drawnow;
    end

    function step(n)
        endFrame = max(startFrame, min(endFrame + n, totalFrames));
        show_frame();
    end

    function jump()
        answer = inputdlg(sprintf('Enter frame number (%d - %d):', startFrame, totalFrames), ...
                          'Jump to Frame', 1, {num2str(endFrame)});
        if ~isempty(answer)
            f = round(str2double(answer{1}));
            if f >= startFrame && f <= totalFrames
                endFrame = f;
                show_frame();
            else
                warndlg(sprintf('Must be between %d and %d.', startFrame, totalFrames), 'Invalid Input');
            end
        end
    end

    function confirm()
        confirmed = true;
        fprintf('End frame confirmed: %d\n', endFrame);
        close(fig);
    end
end
