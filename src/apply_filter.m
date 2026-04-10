function frame = apply_filter(frame, filterType)
% APPLY_FILTER  Apply a pre-processing filter to a video frame.
%
%   filterType options:
%       'none'       - no change (default)
%       'grayscale'  - convert to grayscale (kept as RGB for compatibility)
%       'gaussian'   - Gaussian blur (5x5, sigma 1.5) for noise reduction
%       'sharpen'    - sharpen edges of fast-moving objects

    switch lower(filterType)
        case 'grayscale'
            gray  = rgb2gray(frame);
            frame = cat(3, gray, gray, gray);
        case 'gaussian'
            frame = imfilter(frame, fspecial('gaussian', [5 5], 1.5));
        case 'sharpen'
            frame = imsharpen(frame);
        case 'none'
            % no change
        otherwise
            warning('Unknown filter: "%s" — skipping', filterType);
    end
end
