function export_frames(v, startFrame, endFrame, outputDir, filterType)
% EXPORT_FRAMES  Save frames from startFrame to endFrame as PNGs.
%
%   Inputs:
%       v           - VideoReader object
%       startFrame  - first frame to export (1-indexed)
%       endFrame    - last frame to export
%       outputDir   - path to output directory
%       filterType  - 'none' | 'grayscale' | 'gaussian' | 'sharpen'

    if nargin < 5, filterType = 'none'; end

    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    fprintf('Exporting frames %d to %d with filter: %s\n', startFrame, endFrame, filterType);

    for i = startFrame:endFrame
        v.CurrentTime = (i - 1) / v.FrameRate;
        if hasFrame(v)
            frame = readFrame(v);
            frame = apply_filter(frame, filterType);
            imwrite(frame, fullfile(outputDir, sprintf('frame_%05d.png', i)));
        end
    end

    fprintf('Done! %d frames saved to: %s\n', endFrame - startFrame + 1, outputDir);
end
