function T = survey_capture_integrity(varargin)
%SURVEY_CAPTURE_INTEGRITY  Screen every Tight/Wide raw clip for dropped
%   frames, using a physics bound: during visible motion the top marker
%   cannot move faster than free fall allows. READ-ONLY.
%
%   Per clip (crude red-blob mask, no pipeline needed):
%     nFrames, fps_meta          container facts
%     nMarkerFrames              frames with >=1 marker visible
%     nMovingFrames              frames where the top marker moved > 1 px
%     maxStep_px                 largest per-frame top-marker displacement
%     ffStep_px                  free-fall bound sqrt(2*g*h)/fps in px
%     stepRatio                  maxStep / ffStep   <- the discriminator
%
%   VERDICTS
%     INTACT?   stepRatio <= 1.6 AND nMovingFrames >= 0.5*expected transit
%     BROKEN    stepRatio >= 3   (frames missing between stored frames)
%     CHECK     anything between, or h = 0 (bound is ~0; test inapplicable,
%               judged on nMovingFrames only)
%
%   USAGE
%       T = survey_capture_integrity;                       % Tight + Wide
%       T = survey_capture_integrity('Models',"Wide Model");
%       sortrows(T, 'stepRatio')                            % best first
%
%   OPTIONS
%       'RawRoot'  default D:\ME_GRANULAB\Test Batches\Batch 5
%       'Models'   default ["Tight Model","Wide Model"]
%       'Save'     write CSV next to RawRoot (default true)

opt.RawRoot = 'D:\ME_GRANULAB\Test Batches\Batch 5';
opt.Models  = ["Tight Model","Wide Model"];
opt.Save    = true;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

MMPP = 0.1079;  G = 9.80;                     % m/s^2
rows = {};
for m = 1:numel(opt.Models)
    fdir = fullfile(opt.RawRoot, char(opt.Models(m)), 'GB', 'full');
    V = dir(fullfile(fdir, '*.avi'));
    fprintf('\n--- %s : %d clips ---\n', opt.Models(m), numel(V));
    for i = 1:numel(V)
        vp = fullfile(V(i).folder, V(i).name);
        tok = regexp(V(i).name, '^(\d+)mm', 'tokens', 'once');
        h_m = str2double(tok{1}) / 1000;
        try
            vr = VideoReader(vp);
        catch
            fprintf('  %-16s UNREADABLE\n', V(i).name); continue
        end
        nF = vr.NumFrames;  fps = vr.FrameRate;
        topx = nan(nF,1);
        for k = 1:nF
            f  = read(vr,k);
            bw = f(:,:,1) > 150 & f(:,:,2) < 100;
            cc = bwconncomp(bw);
            if cc.NumObjects > 0
                s = regionprops(cc, 'Centroid', 'Area');
                s = s([s.Area] >= 6);
                if ~isempty(s)
                    c = vertcat(s.Centroid);
                    topx(k) = max(c(:,1));     % marker furthest from bed
                end
            end
        end
        stp   = abs(diff(topx));  stp = stp(isfinite(stp));
        mv    = stp(stp > 1);                              % moving frames
        maxSt = max([mv; 0]);
        ffSt  = sqrt(2*G*max(h_m, 0.01)) * 1000/MMPP / fps; % px/frame bound
        ratio = maxSt / ffSt;
        expTransit = round((0.045/sqrt(2*G*max(h_m,0.01))) * fps); % ~45mm FOV
        if h_m == 0
            verdict = "CHECK(h=0)";
        elseif ratio >= 3
            verdict = "BROKEN";
        elseif ratio <= 1.6 && numel(mv) >= 0.5*expTransit
            verdict = "INTACT?";
        else
            verdict = "CHECK";
        end
        fprintf('  %-16s %5d fr | ratio %6.1f | moving %3d (expect >~%3d) | %s\n', ...
                V(i).name, nF, ratio, numel(mv), expTransit, verdict);
        rows{end+1} = table(opt.Models(m), string(V(i).name), 1000*h_m, nF, fps, ...
            sum(isfinite(topx)), numel(mv), maxSt, ffSt, ratio, verdict, ...
            'VariableNames', {'model','file','h_mm','nFrames','fps_meta', ...
            'nMarkerFrames','nMovingFrames','maxStep_px','ffStep_px', ...
            'stepRatio','verdict'}); %#ok<AGROW>
    end
end
T = vertcat(rows{:});

fprintf('\n--- summary ---\n');
for v = ["INTACT?","CHECK","CHECK(h=0)","BROKEN"]
    fprintf('  %-10s %3d clips\n', v, sum(T.verdict == v));
end
if opt.Save
    p = fullfile(opt.RawRoot, sprintf('capture_integrity_%s.csv', ...
                 char(datetime('now','Format','yyyyMMdd_HHmmss'))));
    writetable(T, p);
    fprintf('wrote: %s\n', p);
end
fprintf(['\nINTACT? clips still need one manual diag_raw_clip look before trust;\n' ...
         'BROKEN clips are only salvageable for dt-free final depth.\n\n']);
end
