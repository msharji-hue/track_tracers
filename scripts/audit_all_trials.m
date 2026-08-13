function T = audit_all_trials(root, varargin)
% AUDIT_ALL_TRIALS  Read-only audit of EVERY *_tracks.mat under 03_RESULTS.
%
% Runs kd_kinematics on all trials with the DEFAULT bed line, re-runs the
% listed SUSPECT trials with an override bed line, and writes one complete
% table (xlsx + csv) to 03_RESULTS/_batch_logs. Writes NO _kin.mat.
%
%   T = audit_all_trials(root)
%   T = audit_all_trials(root,'suspects',{...},'overrideX',10)
%   T = audit_all_trials(root,'modelH',4.00,'shallowBed',2.50)
%   T = audit_all_trials(root,'model','Default','condition','GB/full', ...
%                        'clickSelect',true,'write',false)
%       -> Default geometry only, dry run, click a curve to flag its trialTag
%
% Columns include: condition, dropHeight_mm, trialTag, bedline, impactFrame,
% stopFrame, preFrames, v0_meas, v0_ff, ratio, d_final_cm, d_override_cm,
% dDelta_cm, nPost, flags, passDepth, needsReview, plus anchor diagnostics
% (rodBedMin/Max, anchorInRange) that show whether impactDistPx is reachable.

    p = inputParser;
    addParameter(p,'suspects',{},@iscell);   % bed-line override test; {} = off
    addParameter(p,'overrideX',10,@isnumeric);
    addParameter(p,'modelH',4.00,@isnumeric);        % model height 40 mm
    addParameter(p,'nearFrac',0.93,@isnumeric);      % near-limit band: >=0.93*modelH
    addParameter(p,'shallowBed',2.50,@isnumeric);   % shallow bed depth 25 mm
    addParameter(p,'figures','show',@(x)ischar(x)||isstring(x)); % none|show|save
    addParameter(p,'model','',@(x)ischar(x)||isstring(x));  % ''|Default|Tight|Wide
    addParameter(p,'condition','',@(x)ischar(x)||isstring(x)||iscell(x)); % e.g. 'GB/full'
    addParameter(p,'clickSelect',false,@islogical); % click a curve to flag it
    addParameter(p,'write',true,@islogical);        % false = dry run, no files
    parse(p,varargin{:}); o = p.Results;

    here = fileparts(mfilename('fullpath'));
    addpath(here, fullfile(fileparts(here),'src'));
    calib = get_calibration();  g = calib.g_cm_s2;   % default (reporting/anchor range)
    calibO = get_calibration(o.overrideX);           % optional bed-line override test
    susp = string(o.suspects);

    D = dir(fullfile(root,'03_RESULTS','**','*_tracks.mat'));
    if isempty(D), error('No *_tracks.mat under %s/03_RESULTS', root); end
    n = numel(D);
    fprintf('Auditing %d trials (default bed x=%g; override x=%g on %d suspects)...\n', ...
        n, calib.bedPoint1(1), o.overrideX, numel(susp));

    ws = warning('off','all'); restore = onCleanup(@() warning(ws)); %#ok<NASGU>

    tag=strings(n,1); cond=tag; bedu=tag; flg=tag;
    hmm=nan(n,1); imf=hmm; stf=hmm; pre=hmm; v0m=hmm; v0f=hmm; rat=hmm;
    dfin=hmm; dovr=hmm; ddel=hmm; npo=hmm; rbmin=hmm; rbmax=hmm; inr=false(n,1);
    astp=hmm; htrue=hmm; v0ft=hmm; ratt=hmm; trig=hmm;
    gmdl=strings(n,1);
    pass=false(n,1); rev=false(n,1);
    SER = cell(n,1);            % trimmed series for the per-condition figures

    for k = 1:n
        tp = fullfile(D(k).folder, D(k).name);
        try
            S = load(tp,'meta','tracks'); m = S.meta; tr = S.tracks;
        catch ME
            tag(k)=string(D(k).name); flg(k)="LOADFAIL:"+string(ME.message); rev(k)=true; continue;
        end
        tag(k)  = string(getf(m,'trialTag',''));
        %  geometry from the tag suffix: Tight/Wide carry one, Default does not.
        %  Each foot has its own measured impact trigger -- see get_calibration.
        tl = lower(tag(k));
        if endsWith(tl,"_tight"),    gm = "Tight";
        elseif endsWith(tl,"_wide"), gm = "Wide";
        else,                        gm = "Default";
        end
        gmdl(k) = gm;
        if ~isempty(o.model) && ~strcmpi(gm, o.model), tag(k)=""; continue; end
        cond(k) = string(getf(m,'material','')) + "/" + string(getf(m,'container',''));
        %  condition filter, e.g. 'GB/full' (accepts a cellstr/string array too)
        if ~isempty(o.condition)
            if ~any(strcmpi(cond(k), string(o.condition))), tag(k)=""; continue; end
        end
        hmm(k)  = getf(m,'dropHeight_mm',NaN);
        fps = getf(m,'fps_true',NaN);
        if ~isfinite(fps) && isfield(tr,'fps'), fps = tr.fps; end
        if ~isfinite(fps) || fps < 1000
            flg(k)="BADFPS"; rev(k)=true; continue;
        end

        isSus = any(tag(k)==susp);
        bedu(k) = "x=" + string(calib.bedPoint1(1));

        X=tr.trackedX; Y=tr.trackedY; DT=1/fps;
        calibT = get_calibration([], getf(m,'container',''), gmdl(k)); % per-model trigger
        try
            evalc('kin = kd_kinematics(X,Y,calibT,DT);');
        catch ME
            flg(k)="KINFAIL:"+string(ME.message); rev(k)=true; continue;
        end

        imf(k)=kin.impact_index; stf(k)=kin.stopFrame; pre(k)=kin.impact_index-1;
        npo(k)=kin.stopFrame-kin.impact_index;
        v0m(k)=kin.v0_cm_s;  v0f(k)=sqrt(2*g*(hmm(k)/10));  rat(k)=v0m(k)/v0f(k);
        %  GB/shallow labels are REVERSED on disk -- see src/true_drop_height.m.
        %  Physics must use the TRUE height; the label stays as an identifier.
        htrue(k) = true_drop_height(hmm(k), cond(k));
        v0ft(k)  = sqrt(2*g*(htrue(k)/10));
        ratt(k)  = v0m(k)/v0ft(k);
        trig(k)  = calibT.impactDistPx;
        dfin(k)=kin.d_final_cm;
        astp(k)=getf2(kin,'a_stop_cm_s2',NaN);
        rb = kin.rodBedDist_px(isfinite(kin.rodBedDist_px));
        if ~isempty(rb)
            rbmin(k)=min(rb); rbmax(k)=max(rb);
            inr(k) = calibT.impactDistPx>=rbmin(k) && calibT.impactDistPx<=rbmax(k);
        end

        if isSus
            bedu(k) = "x=" + string(o.overrideX) + " (override)";
            try
                evalc('kinO = kd_kinematics(X,Y,calibO,DT);');
                dovr(k)=kinO.d_final_cm; ddel(k)=dovr(k)-dfin(k);
            catch
                dovr(k)=NaN; ddel(k)=NaN;
            end
        end

        % ── flags: genuine defects only ───────────────────────────────────
        f = strings(0,1);
        dd = dfin(k);
        if ~isfinite(dd),                 f(end+1)="NaNDEPTH"; end
        if isfinite(dd) && dd > o.modelH, f(end+1)="OVERDEPTH"; end
        if isfinite(dd) && dd >= o.nearFrac*o.modelH && dd <= o.modelH, f(end+1)="NEARLIMIT"; end
        if isfinite(dd) && dd < 0.2,      f(end+1)="THIN";  end
        if npo(k) < 5,                    f(end+1)="SHORT"; end
        vp = kin.v(kin.impact_index:kin.stopFrame); dv = diff(vp(:));
        if any(dv > 0.15*abs(kin.v0_cm_s)), f(end+1)="GLITCH"; end
        if ~inr(k),                       f(end+1)="ANCHOR_OOR"; end
        %  Shallow bed is 25 mm deep. Penetration beyond that is FLAGGED for
        %  inspection only -- not treated as a defect and not excluded here.
        if contains(lower(cond(k)),"shallow") && isfinite(dd) && dd > o.shallowBed
            f(end+1)="SHALLOW_EXCEEDS_BED";
        end
        if isempty(f), f="ok"; end
        flg(k) = strjoin(f,"|");

        % trimmed series for figures: impact-5ms .. stop+5ms (keeps memory small)
        if ~strcmpi(o.figures,'none')
            pad = round(0.005*fps);
            lo  = max(1, kin.impact_index-pad);
            hi  = min(numel(kin.t_s), kin.stopFrame+pad);
            idx = lo:hi;
            SER{k} = struct('tms',(kin.t_s(idx)-kin.t_s(kin.impact_index))*1e3, ...
                            'z',kin.z(idx),'v',kin.v(idx),'ag',kin.a_plus_g(idx));
        end

        % NEARLIMIT is informational only: 3.7-3.98 cm is physically attainable,
        % so it is labelled but does NOT force a manual review.
        %  SHALLOW_EXCEEDS_BED is informational and deliberately does NOT
        %  affect passDepth or needsReview.
        pass(k) = isfinite(dd) && dd>0 && dd<=o.modelH;
        rev(k)  = ~pass(k) || contains(flg(k),"GLITCH") ...
                  || contains(flg(k),"THIN") || contains(flg(k),"SHORT");
    end

    keep = tag~="";
    T = table(gmdl(keep),cond(keep),hmm(keep),htrue(keep),tag(keep),bedu(keep),trig(keep), ...
              imf(keep),stf(keep),pre(keep), ...
              v0m(keep),v0f(keep),rat(keep),v0ft(keep),ratt(keep),dfin(keep),astp(keep), ...
              dovr(keep),ddel(keep),npo(keep), ...
              rbmin(keep),rbmax(keep),inr(keep),flg(keep),pass(keep),rev(keep), ...
        'VariableNames',{'model','condition','dropHeight_mm','dropHeight_true_mm','trialTag', ...
        'bedline','impactDistPx_used','impactFrame','stopFrame','preFrames', ...
        'v0_meas','v0_ff','ratio','v0_ff_true','ratio_true','d_final_cm','a_stop_cm_s2', ...
        'd_override_cm','dDelta_cm','nPost','rodBedMin','rodBedMax','anchorInRange', ...
        'flags','passDepth','needsReview'});
    T = sortrows(T,{'condition','dropHeight_mm','trialTag'});

    % ── report ────────────────────────────────────────────────────────────
    fprintf('\n================ ALL-TRIAL AUDIT ================\n');
    fprintf('trials audited      : %d\n', height(T));
    fprintf('pass depth sanity   : %d  (0 < d <= %.2f cm)\n', sum(T.passDepth), o.modelH);
    fprintf('need manual review  : %d\n', sum(T.needsReview));
    fprintf('anchor out of range : %d   <-- impactDistPx=%g unreachable in these trials\n', ...
        sum(~T.anchorInRange), calib.impactDistPx);
    fprintf('\nflag counts:\n');
    fl = ["NaNDEPTH","OVERDEPTH","NEARLIMIT","THIN","SHORT","GLITCH","ANCHOR_OOR","SHALLOW_EXCEEDS_BED"];
    for f = fl, fprintf('  %-12s %d\n', f, sum(contains(T.flags,f))); end

    ov = T(~isnan(T.d_override_cm), :);
    if ~isempty(ov)
        fprintf('\nBED-LINE OVERRIDE (x=%g) on %d suspect trials:\n', o.overrideX, height(ov));
        fprintf('  max |depth change| = %.6f cm   (mean %.6f cm)\n', ...
            max(abs(ov.dDelta_cm)), mean(abs(ov.dDelta_cm),'omitnan'));
        if max(abs(ov.dDelta_cm)) < 1e-9
            fprintf('  => IDENTICAL. A vertical bed-line shift cancels in depth\n');
            fprintf('     (depth is measured relative to the impact frame).\n');
            fprintf('     The real lever is impactDistPx / the impact anchor.\n');
        end
    end

    fprintf('\nv0 ratio (LABELLED height): median %.2f  max %.2f\n', ...
        median(T.ratio,'omitnan'), max(T.ratio));
    fprintf('v0 ratio (TRUE height)    : median %.2f  max %.2f   <- should be <= ~1.1\n', ...
        median(T.ratio_true,'omitnan'), max(T.ratio_true));
    nrev = sum(T.dropHeight_mm ~= T.dropHeight_true_mm);
    fprintf('height labels corrected   : %d trial(s) (GB/shallow)\n', nrev);
    fprintf('\nfigures: one per condition x geometry\n');
    disp(groupsummary(T,{'condition','model'}))

    fprintf('impact trigger used, per geometry:\n');
    disp(groupsummary(T,{'model','impactDistPx_used'}))

    fprintf('\nper-condition pass rate:\n');
    cc = unique(T.condition);
    for i=1:numel(cc)
        s = T.condition==cc(i);
        fprintf('  %-16s %3d trials, %3d pass (%.0f%%), median d=%.2f cm\n', cc(i), sum(s), ...
            sum(T.passDepth(s)), 100*sum(T.passDepth(s))/sum(s), median(T.d_final_cm(s),'omitnan'));
    end

    % ── export ────────────────────────────────────────────────────────────
    logDir = fullfile(root,'03_RESULTS','_batch_logs');
    if ~exist(logDir,'dir'), mkdir(logDir); end
    stamp = datestr(now,'yyyymmdd_HHMMSS');

    % ── per-condition kinematics figures (bad trials highlighted red) ─────
    if o.clickSelect
        setappdata(0,'selectedExcludeTags',{});
        assignin('base','selectedExcludeTags',{});
    end
    if ~strcmpi(o.figures,'none')
        kk = find(keep);
        %  One figure per CONDITION x GEOMETRY. Grouping by condition alone
        %  would overlay Default/Tight/Wide in a single axes.
        grpKey = cond + " | " + gmdl;
        gg = unique(grpKey(kk),'stable');
        for i = 1:numel(gg)
            sel = kk(grpKey(kk)==gg(i) & ~cellfun(@isempty,SER(kk)));
            if isempty(sel), continue; end
            hs = unique(hmm(sel)); hs = hs(isfinite(hs));
            cmapH = parula(max(2,numel(hs)));
            fh = figure('Color','w','Name',['audit — ' char(gg(i))]);
            nBad = 0;
            for j = sel(:)'
                s = SER{j};
                isBad = ~pass(j) || contains(flg(j),"GLITCH");
                if isBad
                    col3 = [0.85 0.1 0.1]; lw = 1.6; nBad = nBad + 1;
                else
                    ci = find(hs==hmm(j),1); if isempty(ci), ci = 1; end
                    col3 = cmapH(ci,:); lw = 0.5;
                end
                nm = char(tag(j));
                subplot(3,1,1); hold on;
                L1 = plot(s.tms, s.z,  '-','Color',col3,'LineWidth',lw,'DisplayName',nm);
                subplot(3,1,2); hold on;
                L2 = plot(s.tms, s.v,  '-','Color',col3,'LineWidth',lw,'DisplayName',nm);
                subplot(3,1,3); hold on;
                L3 = plot(s.z,   s.ag, '-','Color',col3,'LineWidth',lw,'DisplayName',nm);
                if o.clickSelect
                    for L = [L1 L2 L3]
                        L.UserData = struct('origColor',col3,'origWidth',lw);
                        L.ButtonDownFcn = @onPickCurve;
                    end
                end
            end
            subplot(3,1,1); grid on; ylabel('z (cm)'); xlabel('t (ms, 0=impact)');
                yline(o.modelH,'--k',sprintf('model %.2f cm',o.modelH));
                title(sprintf('%s   (%d trials, %d flagged red)', char(gg(i)), numel(sel), nBad), ...
                      'Interpreter','none');
            subplot(3,1,2); grid on; ylabel('v (cm/s)'); xlabel('t (ms, 0=impact)');
            subplot(3,1,3); grid on; ylabel('a+g (cm/s^2)'); xlabel('depth z (cm)');
                xline(o.modelH,'--k');
            if o.clickSelect
                title(subplot(3,1,1), sprintf('%s   (%d trials, %d flagged red)  CLICK a curve to flag it', ...
                      char(gg(i)), numel(sel), nBad),'Interpreter','none');
            end
            if strcmpi(o.figures,'save')
                cname = regexprep(char(gg(i)),'[^\w]','_');
                saveas(fh, fullfile(logDir, sprintf('audit_%s_%s.png', cname, stamp)));
                close(fh);
            end
        end
    end

    if o.write
        csvP  = fullfile(logDir, sprintf('audit_all_trials_%s.csv',  stamp));
        xlsP  = fullfile(logDir, sprintf('audit_all_trials_%s.xlsx', stamp));
        writetable(T, csvP);
        try
            writetable(T, xlsP, 'Sheet','all_trials');
            writetable(T(T.needsReview,:), xlsP, 'Sheet','needs_review');
            fprintf('\nExcel : %s\n', xlsP);
        catch ME
            fprintf('\n(xlsx export failed: %s — CSV still written)\n', ME.message);
        end
        fprintf('CSV   : %s\n', csvP);
    else
        fprintf('\nDRY RUN: no audit tables written.\n');
    end
    if o.clickSelect
        fprintf('CLICK a curve to flag it -> selectedExcludeTags (base ws)\n');
        fprintf('   click again to deselect;  run  show_selected  to print the list\n');
    end

    fprintf('\n--- TRIALS TO INSPECT / EXCLUDE (%d) ---\n', sum(T.needsReview));
    disp(T(T.needsReview, {'condition','dropHeight_mm','trialTag','d_final_cm','nPost','flags'}));
    fprintf('Audit is read-only — no _kin.mat written.\n');
end

% ── helpers ───────────────────────────────────────────────────────────────
function onPickCurve(src, ~)
% Click a plotted curve to flag its trial for manual exclusion.
    nm = src.DisplayName;
    if isempty(nm), return; end
    L = getappdata(0,'selectedExcludeTags'); if isempty(L), L = {}; end
    sib = findall(ancestor(src,'figure'),'Type','line','DisplayName',nm);
    if any(strcmp(L,nm))
        L(strcmp(L,nm)) = [];
        for h = sib(:)'
            if isstruct(h.UserData)
                h.Color = h.UserData.origColor; h.LineWidth = h.UserData.origWidth;
            end
        end
        fprintf('DESELECTED          : %s   (%d selected)\n', nm, numel(L));
    else
        L{end+1} = nm;
        for h = sib(:)'
            h.Color = [0.85 0 0]; h.LineWidth = 2.5; uistack(h,'top');
        end
        fprintf('SELECTED FOR REVIEW : %s   (%d selected)\n', nm, numel(L));
    end
    setappdata(0,'selectedExcludeTags',L);
    assignin('base','selectedExcludeTags',L);
end

function v = getf2(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end
end

function s = default_suspects() %#ok<DEFNU>
%  RETIRED: the bed-line override was shown to be a no-op for depth (a vertical
%  shift cancels in z - z(impact)). Kept only for reference; 'suspects' now
%  defaults to {}.
    s = { ...
    '65mm_T01_as_poured','65mm_T02_as_poured','65mm_T03_as_poured', ...
    '65mm_T04_as_poured','65mm_T05_as_poured','125mm_T01_as_poured', ...
    '125mm_T04_as_poured','165mm_T01_as_poured','165mm_T04_as_poured', ...
    '165mm_T05_as_poured','365mm_T03_as_poured','365mm_T05_as_poured', ...
    '325mm_T03_full','325mm_T04_full','345mm_T02_full','345mm_T03_full','345mm_T04_full', ...
    '25mm_T04_full','125mm_T03_full','185mm_T03_full','205mm_T05_full', ...
    '225mm_T04_full','285mm_T04_full','325mm_T05_full', ...
    '25mm_T03_shallow','25mm_T04_shallow','65mm_T02_shallow','65mm_T03_shallow', ...
    '65mm_T05_shallow','125mm_T03_shallow','125mm_T05_shallow'};
end
