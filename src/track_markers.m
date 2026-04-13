function [trackedX, trackedY] = track_markers(detections, firstFrameCenters, tolerancePx)
% TRACK_MARKERS  Frame-to-frame nearest-neighbour marker assignment.
%
%   Inputs:
%       detections        - cell array of Nx2 center coordinates per frame
%       firstFrameCenters - Mx2 initial marker positions
%       tolerancePx       - max allowed displacement per frame (default 10)
%
%   Outputs:
%       trackedX - M x nFrames matrix of x positions (NaN = not detected)
%       trackedY - M x nFrames matrix of y positions (NaN = not detected)

    if nargin < 3 || isempty(tolerancePx)
        % auto-estimate: use median nearest-neighbour distance in frame 1
        % as a baseline, then allow 2x that as tolerance
        if size(detections{1}, 1) > 1
            D = pdist(detections{1}(:,1:2));
            tolerancePx = 0.5 * min(D);   % half the min inter-marker distance
        else
            tolerancePx = 20;
        end
        fprintf('Auto tolerance: %.1f px\n', tolerancePx);
    end

    % Sort markers by x descending so toe (rightmost) is always marker 1
    [~, sortIdx]      = sort(firstFrameCenters(:,1), 'descend');
    firstFrameCenters = firstFrameCenters(sortIdx, :);

    nFrames  = numel(detections);
    nMarkers = size(firstFrameCenters, 1);
    trackedX = nan(nMarkers, nFrames);
    trackedY = nan(nMarkers, nFrames);
    trackedX(:,1) = firstFrameCenters(:,1);
    trackedY(:,1) = firstFrameCenters(:,2);

    for f = 2:nFrames
        curr = detections{f};
        if isempty(curr), continue; end
        used = false(size(curr,1), 1);

        for m = 1:nMarkers
            last = find(isfinite(trackedX(m,1:f-1)), 1, 'last');
            if isempty(last), continue; end

            dists       = sqrt((curr(:,1)-trackedX(m,last)).^2 + (curr(:,2)-trackedY(m,last)).^2);
            dists(used) = inf;
            [minD, idx] = min(dists);

            if minD <= tolerancePx
                trackedX(m,f) = curr(idx,1);
                trackedY(m,f) = curr(idx,2);
                used(idx)     = true;
            end
        end
    end
end
