% STEP 4 — geometry effects at fixed n = 2/3. No exponent hunting, no
% depth-dependent area. Question: does fixed-2/3 need per-geometry
% coefficients at all, and if so do they follow the hull-area candidate
% normalization? M_hull (d = d0 + b*A_hull^(1/3)*v0^(2/3)) is a candidate
% hypothesis in ORIGINAL depth units, NOT a strict Uehara/KD prediction
% (equal-mass feet violate the sphere law's density assumptions). All
% model comparison happens in cm via AIC/BIC/F; the d/A_hull^(1/3) panel
% is a visualization only.
T = readtable('D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports\master_trials_20260822_215312.csv');
K = T(T.keep_reviewed & ~T.isZeroDrop, :);
K.model = string(K.model);
flag = K.d_final_cm < 1 | K.v0_cm_s < 55;
x = K.v0_cm_s.^(2/3);  y = K.d_final_cm;  N = numel(y);
mdl  = ["Tight","Default","Wide"];
Gm   = double(K.model == mdl);               % N x 3 indicator, column order = mdl
gidx = Gm * [1;2;3];                         % 1/2/3 per trial
Ahull = [2.607 3.495 4.052];                 % locked, order = mdl
h13  = (Ahull(gidx).').^(1/3);               % A_hull^(1/3) per trial (verify Nx1)
if ~iscolumn(h13), h13 = h13(:); end

% ── 0. actual pooled within-(model,height) SD of d ───────────────────
[gid, ~, ~] = findgroups(K.model, K.dropHeight_mm);
ng   = splitapply(@numel, y, gid);
sdg  = splitapply(@std,   y, gid);
sdw  = sqrt( sum((ng-1).*sdg.^2) / sum(ng-1) );   % pooled within-group SD
fprintf('pooled within-group SD of d_final = %.4f cm  (%d groups)\n', sdw, numel(ng));

% ── 1. coefficient ladder, all in cm ─────────────────────────────────
% M0 single | M_hull constrained | M1 a_g | M2 d0_g | M3 full
X    = {[ones(N,1) x], [ones(N,1) h13.*x], [ones(N,1) Gm.*x], [Gm x], [Gm Gm.*x]};
name = ["M0 single","M_hull","M1 a_g","M2 d0_g","M3 full"];
kpar = [2 2 4 4 6];
C = cell(1,5); rss = zeros(1,5);
for m = 1:5
    C{m} = X{m} \ y;  rss(m) = sum((y - X{m}*C{m}).^2);
    fprintf('%-9s k=%d  RMSE=%.4f  AIC=%8.1f  BIC=%8.1f\n', name(m), kpar(m), ...
        sqrt(rss(m)/N), N*log(rss(m)/N)+2*kpar(m), N*log(rss(m)/N)+kpar(m)*log(N));
end
Fst = @(m0,m1) ((rss(m0)-rss(m1))/(kpar(m1)-kpar(m0))) / (rss(m1)/(N-kpar(m1)));
pv  = @(F,d1,d2) 1 - fcdf(F, d1, d2);
fprintf(['F: M0->M1 F=%.2f p=%.4g | M_hull->M1 F=%.2f p=%.4g | ' ...
         'M0->M2 F=%.2f p=%.4g | M1->M3 F=%.2f p=%.4g | M0->M3 F=%.2f p=%.4g\n'], ...
    Fst(1,3), pv(Fst(1,3),2,N-4), Fst(2,3), pv(Fst(2,3),2,N-4), ...
    Fst(1,4), pv(Fst(1,4),2,N-4), Fst(3,5), pv(Fst(3,5),2,N-6), Fst(1,5), pv(Fst(1,5),4,N-6));
fprintf('M_hull vs M0 (both k=2, non-nested): dAIC = %+.1f, dBIC = %+.1f (negative favors M_hull)\n', ...
    (N*log(rss(2)/N)+4) - (N*log(rss(1)/N)+4), ...
    (N*log(rss(2)/N)+2*log(N)) - (N*log(rss(1)/N)+2*log(N)));
% per-geometry mean residuals of M_hull: where does the constraint strain?
rh = y - X{2}*C{2};
offh = arrayfun(@(m) mean(rh(gidx==m)), 1:3);
fprintf('M_hull: d0=%.3f b=%.4f | mean residual T/D/W = %+.4f %+.4f %+.4f cm\n', C{2}, offh);

% ── 2. M1 ratios + matched-v0 differences, one cluster bootstrap ────
c1 = C{3};                                   % [d0 aT aD aW]
fprintf('\nM1: d0=%.3f  a T/D/W = %.4f / %.4f / %.4f\n', c1);
fprintf('point ratios a_D/a_T=%.3f a_W/a_T=%.3f | hypotheses: no-effect 1,1 ; hull 1.103,1.158\n', ...
    c1(3)/c1(2), c1(4)/c1(2));
rng(1); B = 2000;
ug = unique(gid);
vref = [150 200 250];
ratB = nan(B,2);  dDT = nan(B,3);  dWT = nan(B,3);
for b = 1:B
    gs = ug(randi(numel(ug), numel(ug), 1));
    i  = cell2mat(arrayfun(@(u) find(gid==u), gs, 'uni', 0));
    cb = [ones(numel(i),1) Gm(i,:).*x(i)] \ y(i);
    ratB(b,:) = [cb(3)/cb(2), cb(4)/cb(2)];
    for q = 1:3
        dDT(b,q) = (cb(3)-cb(2)) * vref(q)^(2/3);
        dWT(b,q) = (cb(4)-cb(2)) * vref(q)^(2/3);
    end
end
ci = @(v) prctile(v, [2.5 97.5]);
cD = ci(ratB(:,1)); cW = ci(ratB(:,2));
fprintf('cluster CIs: a_D/a_T [%.3f %.3f]  a_W/a_T [%.3f %.3f]\n', cD, cW);
for q = 1:3
    pD = (c1(3)-c1(2))*vref(q)^(2/3);  pW = (c1(4)-c1(2))*vref(q)^(2/3);
    cq = ci(dDT(:,q));  cw = ci(dWT(:,q));
    fprintf('v0=%3d: dD-dT = %+.3f cm [%+.3f %+.3f] (%.2f SD) | dW-dT = %+.3f cm [%+.3f %+.3f] (%.2f SD)\n', ...
        vref(q), pD, cq, pD/sdw, pW, cw, pW/sdw);
end

% ── 3. session check: Default 325/345/365, conditioned on v0 ────────
res = y - X{3}*C{3};
D = K.model=="Default" & ismember(K.dropHeight_mm, [325 345 365]);
day = extractBefore(string(K.capture_datetime), 11);
[sg, sh, sd] = findgroups(K.dropHeight_mm(D), day(D));
disp(table(sh, sd, splitapply(@numel, res(D), sg), splitapply(@mean, res(D), sg), ...
    splitapply(@std, res(D), sg), 'VariableNames', {'h_mm','day','n','res_mu','res_sd'}))

% ── 4. sensitivity without the 10 flagged trials ────────────────────
i = ~flag;  Ni = sum(i);  rs = zeros(1,5);
for m = 1:5
    cs = X{m}(i,:) \ y(i);  rs(m) = sum((y(i) - X{m}(i,:)*cs).^2);
end
bicS = Ni*log(rs/Ni) + kpar*log(Ni);
[~, best] = min(bicS);
c1s = [ones(Ni,1) Gm(i,:).*x(i)] \ y(i);
fprintf('\nsensitivity (%d flagged removed): BIC order best = %s | ratios a_D/a_T=%.3f a_W/a_T=%.3f\n', ...
    sum(flag), name(best), c1s(3)/c1s(2), c1s(4)/c1s(2));

% ── 5. working figure: means±SD (both axes) + M1 curves; raw & normalized
% Right panel is a VISUALIZATION of the hull hypothesis, not a metric.
fig = figure('Position', [80 80 1200 460]);
vv = linspace(60, 290, 200);
for p = 1:2
    subplot(1, 2, p); hold on; grid on
    for m = 1:3
        r = (gidx == m);
        if p == 1, s = 1; else, s = Ahull(m)^(1/3); end
        gi = findgroups(K.dropHeight_mm(r));
        vm = splitapply(@mean, K.v0_cm_s(r), gi);
        vs = splitapply(@std,  K.v0_cm_s(r), gi);
        dm = splitapply(@mean, y(r)./s,      gi);
        ds = splitapply(@std,  y(r)./s,      gi);
        errorbar(vm, dm, ds, ds, vs, vs, 'o');   % yneg,ypos,xneg,xpos
        plot(vv, (c1(1) + c1(1+m).*vv.^(2/3))./s, '-');
    end
    xlabel('v_0 (cm/s)');
    if p == 1
        ylabel('d_{final} (cm)');
        title('raw');
    else
        ylabel('d_{final} / A_{hull}^{1/3}');
        title('hull-normalized (visualization only)');
    end
    legend(["Tight","","Default","","Wide",""], 'Location', 'northwest');
end
exportgraphics(fig, 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures\step4_geometry_working.png', 'Resolution', 150);
