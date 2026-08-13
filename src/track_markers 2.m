function [trackedX, trackedY] = track_markers(detections, firstFrameCenters, tolerancePx)
% TRACK_MARKERS  Frame-to-frame nearest-neighbour marker assignment.
%
%   Inputs:
%       detections        - cell array of Nx2 center coordinates per frame
%       firstFrameCenters - Mx2 initial marker positions
%       tolerancePx       - max match distance from a marker's LAST KNOWN
%                           position. MUST be < one inter-marker spacing, else
%                           a stale reference can match a NEIGHBOUR's detection
%                           (measured spacing ~34.6 px => use ~25 px).
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

    % Sort markers by x descending. NOTE: markers sit on the RIGID ROD; none of
    % them is the foot/toe. With the bed line at LOW x, marker 1 (largest x) is
    % the marker FARTHEST from the bed and the LAST marker is NEAREST the bed.
    % Do not assume marker 1 is 'the toe'.
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
