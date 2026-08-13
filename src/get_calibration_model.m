function calib = get_calibration_model(model, dropHeight_mm, container)
%GET_CALIBRATION_MODEL  Per-model calibration for the Batch 5 foot models.
%
%   calib = get_calibration_model('Default Model',   0, 'full')
%   calib = get_calibration_model('Tight Model',   365, 'full')
%   calib = get_calibration_model('Wide Model',      0, 'full')
%
%   Wraps get_calibration() and overrides only the bed line and the impact
%   trigger, both of which were measured per model before testing. Everything
%   else -- mmPerPx, trackTolerancePx, g, postCapMs -- is inherited unchanged,
%   so this cannot silently alter unrelated settings.
%
%   SIGN. impactDistPx is the SIGNED distance (bedX - x) along the bed normal
%   at contact. The rod marker sits far to the +x side of the line at x = 4, so
%   the value is NEGATIVE. Magnitudes supplied as 370.001 / 376.001 / 409 are
%   applied as -370.001 / -376.001 / -409. A positive value would place the
%   trigger on the side the marker never reaches, so it would never fire and
%   every trial would anchor to a trace endpoint.
%
%   CALIBRATION TABLE (as measured)
%     Default Model   bed (4,0)-(4,32)   trigger -370.001    [d0 trials only]
%     Tight Model     bed (4,0)-(4,24)   trigger -376.001
%                     bed (4,0)-(4,32)   for h = 305/325/345/365 mm
%     Wide Model      bed (4,0)-(4,24)   trigger -409
%                     bed (4,0)-(4,20)   for the d0 (h = 0) trials
%
%   NOTE ON THE BED POINTS. Both points share x = 4, so every variant defines
%   the SAME vertical line; only the segment endpoints differ. The normalised
%   distance used downstream is therefore identical, and lineA/B/C come out
%   proportional. The bed-point variants are recorded for provenance and do
%   not change the computed depth. The trigger is what actually differs
%   between models.

if nargin < 2, dropHeight_mm = NaN; end
if nargin < 3, container = 'full'; end

calib = get_calibration([], container);      % inherit everything else
model = string(model);

% Normalise: strip the word "model" and any spacing, case-insensitively, so
% "Default Model", "default", "Default_Model" and "DEFAULT" all resolve.
key = lower(strtrim(string(model)));
key = erase(key, ["model", "_", " "]);
switch key
    case "default"
        bp2y = 32;
        idp  = -370.001;
    case "tight"
        if ismember(round(dropHeight_mm), [305 325 345 365])
            bp2y = 32;                       % long-drop bed points
        else
            bp2y = 24;
        end
        idp = -376.001;
    case "wide"
        if round(dropHeight_mm) == 0
            bp2y = 20;                       % d0 bed points
        else
            bp2y = 24;
        end
        idp = -409;
    otherwise
        error('get_calibration_model:unknownModel', ...
              'Unknown model "%s" (normalised to "%s"); expected Default, Tight or Wide.', ...
              model, key);
end

calib.model        = char(model);
calib.bedPoint1    = [calib.bedX, 0];
calib.bedPoint2    = [calib.bedX, bp2y];
calib.impactDistPx = idp;

% recompute the implicit line for the recorded endpoints
calib.lineA = calib.bedPoint1(2) - calib.bedPoint2(2);
calib.lineB = calib.bedPoint2(1) - calib.bedPoint1(1);
calib.lineC = calib.bedPoint1(1)*calib.bedPoint2(2) - ...
              calib.bedPoint2(1)*calib.bedPoint1(2);
end