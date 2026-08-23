% STEP 3 — depth-scaling fits on trial-level measured v0. No area terms yet.
T = readtable('D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports\master_trials_20260822_215312.csv');
K = T(T.keep_reviewed & ~T.isZeroDrop, :);
g2 = 2*980;                                        % 2g, cm/s^2
flag = K.d_final_cm < 1 | K.v0_cm_s < 55;          % influence set (kept in fits)
fprintf('flagged for influence checks: %d trials\n', sum(flag));

mdls = ["Tight","Default","Wide","Pooled"];
R = struct();  B = 2000;  rng(1);
for m = 1:4
    if m < 4, r = K.model == mdls(m); else, r = true(height(K),1); end
    x = K.v0_cm_s(r);  y = K.d_final_cm(r);  N = numel(x);
    tag = K.trialTag(r);

    % ── forms ─────────────────────────────────────────────────────────
    [c23, rss23] = lgn(x, y, 2/3);                 % d = d0 + a*v0^(2/3)
    [nf, cf, rssf] = freen(x, y);                  % d = d0 + a*v0^n, n free
    [cl,  rssl ] = lgn(x, y, 1);                   % d = d0 + c*v0
    bp = (x.^(2/3)) \ y;                           % pure 2/3, no intercept
    rssp = sum((y - bp*x.^(2/3)).^2);
    d0mult = bp^(3/2) * sqrt(g2);                  % multiplicative-form d0: d=(d0^2 H)^(1/3)

    % ── bootstrap n (trials, and clusters by (model,height)) ─────────
    nb = nan(B,1); nbc = nan(B,1);
    [gid,~,~] = findgroups(K.model(r), K.dropHeight_mm(r)); ug = unique(gid);
    for b = 1:B
        i = randi(N,N,1);                nb(b)  = freen(x(i), y(i));
        gs = ug(randi(numel(ug),numel(ug),1));
        i = cell2mat(arrayfun(@(u) find(gid==u), gs, 'uni',0));
        nbc(b) = freen(x(i), y(i));
    end
    ci  = prctile(nb, [2.5 97.5]);  cic = prctile(nbc, [2.5 97.5]);

    % ── group-mean weighted fit (literature-style cross-check) ───────
    v0m = splitapply(@mean,x,gid); dm = splitapply(@mean,y,gid); nm = splitapply(@numel,x,gid);
    W = diag(nm); Xg = [ones(size(v0m)) v0m.^(2/3)];
    cg = (Xg'*W*Xg) \ (Xg'*W*dm);
    ng = freen_w(v0m, dm, nm);

    % ── information criteria ─────────────────────────────────────────
    aic = @(rss,k) N*log(rss/N) + 2*k;   bic = @(rss,k) N*log(rss/N) + k*log(N);
    fprintf('\n=== %s (N=%d) ===\n', mdls(m), N);
    fprintf('fixed 2/3 : d0=%6.3f  a=%7.4f   RMSE=%.3f  R2=%.3f  AIC=%7.1f BIC=%7.1f\n', ...
        c23(1), c23(2), sqrt(rss23/N), 1-rss23/sum((y-mean(y)).^2), aic(rss23,2), bic(rss23,2));
    fprintf('free n    : d0=%6.3f  a=%7.4f  n=%.3f  CI_trial [%.3f %.3f]  CI_cluster [%.3f %.3f]  RMSE=%.3f  AIC=%7.1f BIC=%7.1f\n', ...
        cf(1), cf(2), nf, ci, cic, sqrt(rssf/N), aic(rssf,3), bic(rssf,3));
    fprintf('linear    : d0=%6.3f  c=%7.4f   RMSE=%.3f  AIC=%7.1f BIC=%7.1f\n', ...
        cl(1), cl(2), sqrt(rssl/N), aic(rssl,2), bic(rssl,2));
    fprintf('pure 2/3  : b=%.4f -> multiplicative d0=%.3f cm   RMSE=%.3f  AIC=%7.1f BIC=%7.1f\n', ...
        bp, d0mult, sqrt(rssp/N), aic(rssp,2), bic(rssp,2));
    fprintf('group-mean: fixed-2/3 d0=%.3f a=%.4f | free n=%.3f  (weighted by n_g, %d groups)\n', cg(1), cg(2), ng, numel(nm));

    % ── profile RSS(n): the d0-n degeneracy, and data for Fig 3 ──────
    ngrid = 0.2:0.01:1.4;
    rssn = arrayfun(@(nn) second(@() lgn(x,y,nn)), ngrid);
    R.(mdls(m)) = struct('c23',c23,'nf',nf,'cf',cf,'cl',cl,'bp',bp,'d0mult',d0mult, ...
        'ci',ci,'cic',cic,'ngrid',ngrid,'rssn',rssn,'x',x,'y',y,'nb',nb);

    % ── influence: leave-one-out over flagged trials in this set ─────
    fl = find(flag(r));
    for k = fl(:)'
        i = setdiff(1:N, k);
        [nk, ck] = freen(x(i), y(i));
        if abs(nk-nf) > 0.02 || abs(ck(1)-cf(1)) > 0.05
            fprintf('  LOO %-26s : n %.3f->%.3f  d0 %.3f->%.3f\n', tag{k}, nf, nk, cf(1), ck(1));
        end
    end
end

% ── sensitivity: all forms without the flagged set (pooled) ───────────
x = K.v0_cm_s(~flag); y = K.d_final_cm(~flag);
[c23s,~] = lgn(x,y,2/3); [nfs,cfs] = freen(x,y);
fprintf('\nsensitivity (pooled, %d flagged removed): fixed-2/3 d0=%.3f a=%.4f | free n=%.3f d0=%.3f\n', ...
    sum(flag), c23s(1), c23s(2), nfs, cfs(1));

% ── working figures ───────────────────────────────────────────────────
P = R.Pooled;
figure('Position',[80 80 1500 420]);
subplot(1,3,1); gscatter(K.v0_cm_s, K.d_final_cm, K.model); hold on
vv = linspace(40,300,200);
plot(vv, P.c23(1)+P.c23(2)*vv.^(2/3), 'k-', vv, P.cf(1)+P.cf(2)*vv.^P.nf, 'k--');
xlabel('v_0 (cm/s)'); ylabel('d (cm)'); title('fixed 2/3 (solid) vs free n (dashed)'); grid on
subplot(1,3,2);
res = K.d_final_cm - (P.c23(1)+P.c23(2)*K.v0_cm_s.^(2/3));
gscatter(K.v0_cm_s, res, K.model); yline(0); grid on
xlabel('v_0 (cm/s)'); ylabel('residual, fixed 2/3 (cm)'); title('curvature check');
subplot(1,3,3); plot(P.ngrid, P.rssn/min(P.rssn), '-'); hold on
xline(2/3,'--'); xline(P.nf); xline(1,':'); grid on
xlabel('n'); ylabel('RSS(n)/RSS_{min}'); title(sprintf('profile: n=%.3f, CI [%.2f %.2f]', P.nf, P.ci));
exportgraphics(gcf, 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures\step3_fits_working.png', 'Resolution', 200);

% ═══ local functions ══════════════════════════════════════════════════
function [c, rss] = lgn(x, y, n)
    X = [ones(size(x)) x.^n];  c = X\y;  rss = sum((y - X*c).^2);
end
function [n, c, rss] = freen(x, y)
    n = fminbnd(@(nn) second(@() lgn(x,y,nn)), 0.2, 1.4);
    [c, rss] = lgn(x, y, n);
end
function n = freen_w(x, y, w)
    f = @(nn) wrss(x,y,w,nn);  n = fminbnd(f, 0.2, 1.4);
end
function r = wrss(x,y,w,nn)
    X = [ones(size(x)) x.^nn];  W = diag(w);
    c = (X'*W*X)\(X'*W*y);  r = sum(w.*(y-X*c).^2);
end
function v = second(f)
    [~, v] = f();
end
