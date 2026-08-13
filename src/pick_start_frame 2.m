function startFrame = pick_start_frame(v)
% PICK_START_FRAME  Visually browse frames and confirm your start frame.
%   Shows each frame in a figure with buttons:
%       [<<]  step back 10 frames
%       [<]   step back 1 frame
%       [>]   step forward 1 frame
%       [>>]  step forward 10 frames
%       [Jump to Frame]  type a specific frame number
%       [✓ Confirm Start Frame]  lock in your choice

    totalFrames = floor(v.Duration * v.FrameRate);
    startFrame  = 1;
    confirmed   = false;

    % ── Build figure ─────────────────────────────────────────
    fig = figure('Name', 'Pick Start Frame', ...
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
    uicontrol('Style', 'pushbutton', 'String', '✓ Confirm Start Frame', ...
        'Position', [460 20 200 40], ...
        'BackgroundColor', [0.2 0.7 0.2], 'ForegroundColor', 'white', ...
        'FontWeight', 'bold', ...
        'Callback', @(~,~) confirm());

    show_frame();
    waitfor(fig);

    if ~confirmed
        error('No frame selected — figure was closed without confirming.');
    end

    % ── Nested helpers ────────────────────────────────────────
    function show_frame()
        v.CurrentTime = (startFrame - 1) / v.FrameRate;
        imshow(readFrame(v), 'Parent', ax);
        title(ax, sprintf('Frame  %d / %d', startFrame, totalFrames), ...
              'FontSize', 14);
        drawnow;
    end

    function step(n)
        startFrame = max(1, min(startFrame + n, totalFrames));
        show_frame();
    end

    function jump()
        answer = inputdlg(sprintf('Enter frame number (1 - %d):', totalFrames), ...
                          'Jump to Frame', 1, {num2str(startFrame)});
        if ~isempty(answer)
            f = round(str2double(answer{1}));
            if f >= 1 && f <= totalFrames
                startFrame = f;
                show_frame();
            else
                warndlg('Frame number out of range.', 'Invalid Input');
            end
        end
    end

    function confirm()
        confirmed = true;
        fprintf('Start frame confirmed: %d\n', startFrame);
        close(fig);
    end
end
