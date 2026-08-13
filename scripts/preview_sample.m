function T = preview_sample(root, varargin)
% PREVIEW_SAMPLE  Representative pre-batch sweep. Samples 1-2 trials per drop
% height for each condition, computes kinematics (READ-ONLY, writes nothing),
% prints a QC-flag table, and draws one figure per condition (z, v, a+g overlaid
% and colour-coded by drop height). Use before the full Step-7 batch.
%
%   T = preview_sample(root)                      % 1 trial per height
%   T = preview_sample(root,'perHeight',2)        % 2 per height
%   T = preview_sample(root,'figures','none')     % table only, no plots

    p = inputParser;
    addParameter(p,'perHeight',1,@isnumeric);
    addParameter(p,'figures','show',@(x)ischar(x)||isstring(x));
    parse(p,varargin{:}); o = p.Results;

    here = fileparts(mfilename('fullpath'));
    addpath(here, fullfile(fileparts(here),'src'));
    calib = get_calibration(); g = calib.g_cm_s2;

    D = dir(fullfile(root,'03_RESULTS','**','*_tracks.mat'));
    if isempty(D), error('No *_tracks.mat under %s/03_RESULTS', root); end

    % ── light pass: condition / height / trialNum for sampling ────────────
    n = numel(D);
    cond=strings(n,1); h=nan(n,1); tn=nan(n,1); tag=strings(n,1); path=strings(n,1);
    for k=1:n
        tp=fullfile(D(k).folder,D(k).name);
        try s=load(tp,'meta'); m=s.meta; catch, continue; end
        cond(k)=string(getf(m,'material',''))+"/"+string(getf(m,'container',''));
        h(k)=getf(m,'dropHeight_mm',NaN); tn(k)=getf(m,'trialNum',NaN);
        tag(k)=string(getf(m,'trialTag','')); path(k)=string(tp);
    end
    ok = tag~=""; cond=cond(ok); h=h(ok); tn=tn(ok); tag=tag(ok); path=path(ok);

    % ── choose perHeight trials for each (condition, height) ──────────────
    pick = false(numel(tag),1);
    for c = unique(cond)'
        for hh = unique(h(cond==c))'
            idx = find(cond==c & h==hh);
            [~,ord] = sort(tn(idx)); idx = idx(ord);
            pick(idx(1:min(o.perHeight,numel(idx)))) = true;
        end
    end
    sel = find(pick);
    fprintf('Sampling %d trials (%d conditions) ...\n', numel(sel), numel(unique(cond)));

    % ── compute kinematics for the sample (quiet) ─────────────────────────
    warnState = warning('off','kd_kinematics:shortFall');
    restore = onCleanup(@() warning(warnState)); %#ok<NASGU>
    K = cell(numel(sel),1);
    rows = strings(0,1);
    Cn=strings(0,1); Hh=[]; Tg=strings(0,1); Pre=[]; V0=[]; Vff=[]; Rat=[]; Dd=[]; Npost=[]; Flg=strings(0,1);
    for j = 1:numel(sel)
        i = sel(j);
        try
            S=load(path(i),'meta','tracks'); m=S.meta; tr=S.tracks;
            fps=getf(m,'fps_true',NaN); if ~isfinite(fps)&&isfield(tr,'fps'), fps=tr.fps; end
            X=tr.trackedX; Y=tr.trackedY; DT=1/fps;
            evalc('kin = kd_kinematics(X,Y,calib,DT);');
        catch ME
            Cn(end+1,1)=cond(i); Hh(end+1,1)=h(i); Tg(end+1,1)=tag(i);
            Pre(end+1,1)=NaN; V0(end+1,1)=NaN; Vff(end+1,1)=NaN; Rat(end+1,1)=NaN;
            Dd(end+1,1)=NaN; Npost(end+1,1)=NaN; Flg(end+1,1)="ERROR:"+string(ME.message); K{j}=[];
            continue;
        end
        K{j} = struct('kin',kin,'cond',cond(i),'h',h(i),'tag',tag(i));

        pre  = kin.impact_index-1;
        npost= kin.stopFrame - kin.impact_index;
        vff  = sqrt(2*g*(h(i)/10));
        ratio= kin.v0_cm_s/vff;
        vpost= kin.v(kin.impact_index:kin.stopFrame); dv = diff(vpost(:));
        %  NOTE: short pre-impact fall (pre<=1) and v0 differing from sqrt(2gh)
        %  are NOT flags. The camera FOV sits on the bed, so little/no fall is
        %  captured by design; and the rails carry friction, so v0 < sqrt(2gh)
        %  is the expected physics (the reason g_eff exists). Both are reported
        %  as columns (preFrames, ratio) for comparison only.
        f = strings(0,1);
        if kin.d_final_cm<0.2,     f(end+1)="THIN";   end
        if npost<5,                f(end+1)="SHORT";  end
        if any(dv > 0.15*abs(kin.v0_cm_s)), f(end+1)="GLITCH"; end
        if isempty(f), f="ok"; end

        Cn(end+1,1)=cond(i); Hh(end+1,1)=h(i); Tg(end+1,1)=tag(i);
        Pre(end+1,1)=pre; V0(end+1,1)=kin.v0_cm_s; Vff(end+1,1)=vff; Rat(end+1,1)=ratio;
        Dd(end+1,1)=kin.d_final_cm; Npost(end+1,1)=npost; Flg(end+1,1)=strjoin(f,"|");
    end

    T = table(Cn,Hh,Tg,Pre,V0,Vff,Rat,Dd,Npost,Flg, ...
        'VariableNames',{'condition','dropHeight_mm','trialTag','preFrames', ...
        'v0_meas','v0_ff','ratio','d_final_cm','nPost','flags'});
    T = sortrows(T,{'condition','dropHeight_mm'});
    fprintf('\n'); disp(T);

    % ── reported-for-comparison summary (NOT pass/fail) ───────────────────
    fprintf('Pre-impact fall captured: median %d frames (range %d-%d) — set by FOV, informational.\n', ...
        round(median(Pre,'omitnan')), min(Pre), max(Pre));
    r = Rat(isfinite(Rat));
    fprintf('v0_meas / sqrt(2gh): median %.2f  (range %.2f-%.2f)\n', median(r), min(r), max(r));
    fprintf('   ratio < 1 is expected — rail friction reduces the effective drop acceleration.\n');
    fprintf('   implied g_eff/g from the median ratio: %.2f  (g_eff ~ %.0f cm/s^2)\n', ...
        median(r)^2, median(r)^2*g);
    if any(r > 1.05)
        fprintf('   note: %d sampled trial(s) have ratio > 1, which friction alone cannot explain\n', sum(r>1.05));
        fprintf('         (few pre-impact frames make the impact-frame v0 read noisy).\n');
    end

    bad = T(T.flags~="ok", :);
    if isempty(bad)
        fprintf('\nNo defect flags (THIN / SHORT / GLITCH) in the sample.\n');
    else
        fprintf('\nTrials worth a look:\n'); disp(bad);
    end

    % ── one figure per condition: z, v, a+g vs depth, coloured by height ──
    if ~strcmpi(o.figures,'none')
        conds = unique(Cn,'stable');
        for c = conds'
            mem = find(Cn==c & ~cellfun(@isempty,K(1:numel(Cn))));
            if isempty(mem), continue; end
            hs = unique(Hh(mem)); cmap = parula(max(2,numel(hs)));
            f = figure('Color','w','Name',['preview — ' char(c)]);
            leg = strings(0,1);
            for mi = mem'
                kk = K{mi}; if isempty(kk), continue; end
                ci = cmap(find(hs==kk.h,1),:);
                subplot(3,1,1); hold on; plot(kk.kin.t_s*1e3, kk.kin.z,'-','Color',ci);
                subplot(3,1,2); hold on; plot(kk.kin.t_s*1e3, kk.kin.v,'-','Color',ci);
                subplot(3,1,3); hold on; plot(kk.kin.z, kk.kin.a_plus_g,'-','Color',ci);
                leg(end+1,1)=string(kk.h)+" mm";
            end
            subplot(3,1,1); grid on; ylabel('z (cm)'); xlabel('t (ms)');
                title(char(c),'Interpreter','none'); legend(leg,'Location','eastoutside');
            subplot(3,1,2); grid on; ylabel('v (cm/s)'); xlabel('t (ms)');
            subplot(3,1,3); grid on; ylabel('a+g (cm/s^2)'); xlabel('depth z (cm)');
        end
    end
    fprintf('\nPreview is read-only — no _kin.mat written. Full batch will process all trials.\n');
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end
end