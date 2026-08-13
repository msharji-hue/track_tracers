function diag_capture(root)
% DIAG_CAPTURE  Is the impact event actually present in the Tight/Wide tracks?
%
%   Two questions, since raw z for Tight/Wide is a step function rather than a
%   penetration curve:
%     (1) how long is each recording, in frames AND in ms?
%     (2) do all 8 markers survive through the event, or are they lost?
%
%   diag_capture('D:\ME_GRANULAB\JerboaImpact')

    if nargin < 1, root = 'D:\ME_GRANULAB\JerboaImpact'; end
    D = dir(fullfile(root,'03_RESULTS','**','*_tracks.mat'));
    n = numel(D);

    mdl=strings(n,1); tg=mdl; nfr=nan(n,1); fps=nan(n,1); durms=nan(n,1);
    nAll8=nan(n,1); fracAll8=nan(n,1); maxGap=nan(n,1);

    for k = 1:n
        tp = fullfile(D(k).folder, D(k).name);
        try, S = load(tp,'meta','tracks'); catch, continue; end
        t = string(erase(D(k).name,'_tracks.mat'));
        if ~contains(lower(t),'full'), continue; end        % GB/full only
        tg(k)  = t;
        if endsWith(lower(t),"_tight"), mdl(k)="Tight";
        elseif endsWith(lower(t),"_wide"), mdl(k)="Wide";
        else, mdl(k)="Default"; end

        X = S.tracks.trackedX;
        nfr(k) = size(X,2);
        f = resolve_fps(tp, S.meta, S.tracks);
        fps(k) = f;
        durms(k) = nfr(k)/f*1e3;

        good = sum(isfinite(X),1);            % markers visible per frame
        nAll8(k)   = sum(good == size(X,1));
        fracAll8(k)= nAll8(k)/nfr(k);
        % longest run of frames with fewer than all markers
        bad = good < size(X,1);
        d = diff([0 bad 0]);
        starts = find(d==1); ends = find(d==-1)-1;
        if isempty(starts), maxGap(k)=0; else, maxGap(k)=max(ends-starts+1); end
    end

    keep = tg~="";
    T = table(mdl(keep),tg(keep),fps(keep),nfr(keep),durms(keep), ...
              fracAll8(keep),maxGap(keep), ...
        'VariableNames',{'model','trialTag','fps','nFrames','duration_ms', ...
                         'fracAllMarkers','maxDropoutRun'});

    fprintf('\n===== RECORDING LENGTH by model =====\n');
    disp(varfun(@(x)[numel(x) min(x) median(x) max(x)], T, ...
         'InputVariables','nFrames','GroupingVariables','model'))
    fprintf('===== DURATION (ms) by model =====\n');
    disp(varfun(@(x)[min(x) median(x) max(x)], T, ...
         'InputVariables','duration_ms','GroupingVariables','model'))
    fprintf('  Default penetration alone lasts ~30 ms. Any recording shorter\n');
    fprintf('  than that cannot contain a full impact event.\n');

    short = T(T.duration_ms < 40, :);
    fprintf('\ntrials shorter than 40 ms: %d of %d\n', height(short), height(T));
    if ~isempty(short), disp(groupsummary(short,'model')); end

    fprintf('\n===== MARKER COMPLETENESS by model =====\n');
    fprintf('fraction of frames with ALL markers tracked:\n');
    disp(varfun(@(x)[min(x) median(x) max(x)], T, ...
         'InputVariables','fracAllMarkers','GroupingVariables','model'))
    fprintf('longest run of frames with a marker missing:\n');
    disp(varfun(@(x)[median(x) max(x)], T, ...
         'InputVariables','maxDropoutRun','GroupingVariables','model'))

    figure('Color','w','Position',[80 80 900 400]);
    subplot(1,2,1);
    boxchart(categorical(T.model), T.duration_ms); grid on;
    yline(30,'r--','~30 ms event'); ylabel('recording duration (ms)');
    title('Is the event even captured?');
    subplot(1,2,2);
    boxchart(categorical(T.model), T.fracAllMarkers); grid on;
    ylabel('fraction of frames with all markers'); title('Marker completeness');

    assignin('base','CAP',T);
    fprintf('\ntable written to base workspace as CAP\n');
end
