function K = load_kinematics_set(root, varargin)
% LOAD_KINEMATICS_SET  Load saved kinematics into one table for analysis.
%
%   K = load_kinematics_set(root)
%       'model'      ''|'Default'|'Tight'|'Wide'   geometry from tag suffix
%       'condition'  ''|'GB/full'|...              material/container
%       'exclude'    {}                            trialTags to drop
%       'series'     false                         also return z/v/a+g arrays
%
%   Reads <container>/kinematics/*_kin_scalars.csv (and *_kin.mat when
%   'series' is true). Returns one row per trial with, at minimum:
%     trialTag, model, condition, dropHeight_mm, fps, impactFrame, stopFrame,
%     v0_cm_s, d_final_cm, a_stop_cm_s2, impactDistPx, bedX
%
%   ZERO-DROP v0. The protocol defines v0 = 0 at h = 0: a d0 trial is released
%   from contact, so there is no fall and no impact speed. v0_cm_s is therefore
%   set to exactly 0 on those rows. The value kd_kinematics measured is kept as
%   v0_meas_cm_s for QA -- it is the velocity peak of a release transient, and a
%   large one flags a trial that was dropped rather than released. Never fit
%   against v0_meas_cm_s; use it only to check the zero-drop trials.
%
%   dropHeight_mm is the physical drop height and is the column every height
%   axis should use. trialTag is an identifier only.
%
%   Replaces load_default_gb, which hardcoded GB/full + Default geometry and
%   one manual exclusion list. The inclusion rules that were baked into that
%   function are available here as options, so the analysis scripts no longer
%   each carry their own copy:
%
%   OPTIONS
%       'model'       '' (all) | 'Default' | 'Tight' | 'Wide'
%       'condition'   '' (all) | 'GB/full' | 'GB/shallow' | 'CHIN/dense' | ...
%       'exclude'     cellstr/string array of trialTags to drop outright
%       'series'      false; when true also loads kin.z / kin.v / kin.t_s and
%                     returns them as a cell column per trial
%       'dropNaNDepth' true;  drop h > 0 trials with non-finite d_final_cm
%       'dropGlitch'   true;  drop h > 0 trials failing the monotonic-|v| test
%       'depthCut'     [];    drop h > 0 trials deeper than this (cm)
%       'verbose'      true;  print the per-rule breakdown
%
%   GEOMETRY. Two tag conventions coexist: the validated drop trials carry no
%   model suffix (165mm_T03_full) while the newer zero-drop trials do
%   (0mm_T01_full_default). Tight and Wide always carry _tight / _wide. So an
%   unsuffixed tag is Default, not "unknown".
%
%   ZERO-DROP TRIALS (h = 0) are returned with isZeroDrop true and are exempt
%   from the automatic rules, which assume a fall and a deceleration that a
%   released-from-rest trial does not have. Manual exclusions still apply.
%
%   The returned table always carries the diagnostic columns keep and reason,
%   so an exclusion can be traced rather than silently applied. Rows failing a
%   rule are REMOVED from K; pass 'verbose' to see the counts.

opt.model        = '';
opt.condition    = '';
opt.exclude      = strings(0,1);
opt.series       = false;
opt.dropNaNDepth = true;
opt.dropGlitch   = true;
opt.depthCut     = [];
opt.verbose      = true;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

exclude = string(opt.exclude(:));

fprintf('\n=== load_kinematics_set ===\n');

% ── gather candidates ────────────────────────────────────────────────────
F = dir(fullfile(root,'03_RESULTS','**','*_kin_scalars.csv'));
F = F(~[F.isdir]);
if isempty(F)
    error('load_kinematics_set:noFiles', ...
          'No *_kin_scalars.csv under %s', fullfile(root,'03_RESULTS'));
end

n = numel(F);
tag   = strings(n,1);  cond  = strings(n,1);  mdl = strings(n,1);
hLab  = nan(n,1);      d     = nan(n,1);      v0  = nan(n,1);
tstop = nan(n,1);      phi   = nan(n,1);      fps = nan(n,1);
impF  = nan(n,1);      stopF = nan(n,1);      aSt = nan(n,1);
idp   = nan(n,1);      bedX  = nan(n,1);      kp  = strings(n,1);

for i = 1:n
    T = readtable(fullfile(F(i).folder,F(i).name));
    if height(T) < 1, continue; end
    tag(i)   = string(erase(F(i).name,'_kin_scalars.csv'));
    cond(i)  = local_str(T,'condition');
    mdl(i)   = local_model_from_tag(tag(i));
    d(i)     = local_num(T,'d_final_cm');
    v0(i)    = abs(local_num(T,'v0_cm_s'));
    tstop(i) = local_num(T,'t_stop_s');
    phi(i)   = local_num(T,'phi');
    fps(i)   = local_num(T,'fps');
    impF(i)  = local_num(T,'impact_index');
    stopF(i) = local_num(T,'stopFrame');
    aSt(i)   = local_num(T,'a_stop_cm_s2');
    idp(i)   = local_num(T,'impactDistPx');
    bedX(i)  = local_num(T,'bedX');
    hLab(i)  = local_num(T,'dropHeight_mm');
    kp(i)    = string(strrep(fullfile(F(i).folder,F(i).name), ...
                             '_kin_scalars.csv','_kin.mat'));
end

ok    = tag ~= "";
tag=tag(ok); cond=cond(ok); mdl=mdl(ok); hLab=hLab(ok); d=d(ok); v0=v0(ok);
tstop=tstop(ok); phi=phi(ok); fps=fps(ok); impF=impF(ok); stopF=stopF(ok);
aSt=aSt(ok); idp=idp(ok); bedX=bedX(ok); kp=kp(ok);

K = table(tag, mdl, cond, hLab, fps, impF, stopF, ...
          v0, d, tstop, aSt, phi, idp, bedX, kp, ...
    'VariableNames', {'trialTag','model','condition','dropHeight_mm', ...
        'fps','impactFrame','stopFrame', ...
        'v0_cm_s','d_final_cm','t_stop_s','a_stop_cm_s2','phi', ...
        'impactDistPx','bedX','kinPath'});

% ── duplicate trialTags are a hard error ─────────────────────────────────
% One trial must contribute exactly one row. A tag appearing twice means two
% _kin_scalars.csv files describe the same trial -- typically a stale copy left
% under an old batch or model folder after a re-run or a rename. Silently
% keeping both double-weights that trial in every per-height mean and every
% fit, so this stops rather than warns: the caller cannot tell from the results
% that it happened, and the fix is a deliberate deletion, not a filter.
[uTags, ~, uIdx] = unique(K.trialTag);
counts = accumarray(uIdx, 1);
dupTags = uTags(counts > 1);
if ~isempty(dupTags)
    msg = sprintf(['%d trialTag(s) appear more than once. One trial must ' ...
                   'produce exactly one row.\n\n'], numel(dupTags));
    for i = 1:numel(dupTags)
        rows = find(K.trialTag == dupTags(i));
        msg  = [msg sprintf('  %s  (%d copies)\n', dupTags(i), numel(rows))]; %#ok<AGROW>
        for r = rows(:).'
            % kinPath was derived from the scalars CSV; invert that to name the
            % file the user actually has to delete.
            csvPath = strrep(K.kinPath(r), '_kin.mat', '_kin_scalars.csv');
            msg = [msg sprintf('      %s\n', csvPath)]; %#ok<AGROW>
        end
    end
    msg = [msg sprintf(['\nDelete the stale copies (keep exactly one per trial) ' ...
                        'and re-run.\nDo not filter them out here: which copy is ' ...
                        'current is a decision only you can make.\n'])];
    error('load_kinematics_set:duplicateTrialTag', '%s', msg);
end

K.isZeroDrop = K.dropHeight_mm == 0;

% Zero-drop protocol: a d0 trial is released FROM CONTACT, so the impact speed
% is defined to be zero. What kd_kinematics reports for those trials is the
% velocity peak of a release transient, not a fall, so it is a QA number rather
% than a physical impact speed. Keep the measured value under v0_meas_cm_s and
% set the physical column to the protocol value.
K.v0_meas_cm_s = K.v0_cm_s;
K.v0_cm_s(K.isZeroDrop) = 0;

% ── phi: take the CURRENT value, not the one baked into the CSV ──────────
% _kin_scalars.csv carries meta.phi as it stood at TRACKING time, which for
% this dataset is the retired single-measurement figure (0.629 / 0.636 /
% 0.276 / 0.409). get_substrate_properties holds the finalised five-prep
% means (0.624 / 0.643 / 0.280 / 0.402), so re-read it here and overwrite.
K.phi_csv = K.phi;                       % kept for provenance
for i = 1:height(K)
    parts = split(K.condition(i), "/");
    if numel(parts) ~= 2, continue; end
    try
        sub = get_substrate_properties(char(parts(1)), char(parts(2)));
        if isstruct(sub) && isfield(sub,'phi') && isfinite(sub.phi)
            K.phi(i) = sub.phi;
        end
    catch
        % leave the CSV value in place if the lookup fails
    end
end
nPhi = sum(isfinite(K.phi) & isfinite(K.phi_csv) & abs(K.phi-K.phi_csv) > 1e-6);
if nPhi > 0
    fprintf('  phi refreshed from get_substrate_properties : %d row(s)\n', nPhi);
end

fprintf('kin_scalars found         : %d\n', height(K));

% ── model / condition filter ─────────────────────────────────────────────
if ~isempty(opt.model)
    K = K(strcmpi(K.model, string(opt.model)), :);
    fprintf('  model == %-10s      : %d\n', string(opt.model), height(K));
end
if ~isempty(opt.condition)
    K = K(K.condition == string(opt.condition), :);
    fprintf('  condition == %-10s  : %d\n', string(opt.condition), height(K));
end
if isempty(K)
    error('load_kinematics_set:noneSelected', ...
          'No trials left after the model/condition filter.');
end

% ── rules ────────────────────────────────────────────────────────────────
exManual = ismember(K.trialTag, exclude);

exNaN = false(height(K),1);
if opt.dropNaNDepth
    exNaN = ~K.isZeroDrop & ~isfinite(K.d_final_cm);
end

exGlitch = false(height(K),1);
if opt.dropGlitch
    for i = 1:height(K)
        if K.isZeroDrop(i) || exManual(i) || exNaN(i), continue; end
        exGlitch(i) = local_isglitch(K.kinPath(i));
    end
end

exCut = false(height(K),1);
if ~isempty(opt.depthCut)
    exCut = ~K.isZeroDrop & K.d_final_cm > opt.depthCut;
end

K.keep   = ~(exManual | exNaN | exGlitch | exCut);
K.reason = strings(height(K),1);
K.reason(exManual) = "manual";
K.reason(exNaN)    = "NaNDEPTH";
K.reason(exGlitch) = "GLITCH";
K.reason(exCut)    = "depth cut";

if opt.verbose
    fprintf('\n--- exclusions applied ---\n');
    fprintf('  manual   : %2d of %d listed  (%d listed tags matched no file)\n', ...
        sum(exManual), numel(exclude), sum(~ismember(exclude, K.trialTag)));
    miss = exclude(~ismember(exclude, K.trialTag));
    if ~isempty(miss), fprintf('             %s\n', miss); end
    if opt.dropNaNDepth, fprintf('  NaNDEPTH : %2d  (h > 0 only)\n', sum(exNaN));    end
    if opt.dropGlitch,   fprintf('  GLITCH   : %2d  (h > 0 only)\n', sum(exGlitch)); end
    if ~isempty(opt.depthCut)
        fprintf('  depth cut: %2d  (> %.2f cm)\n', sum(exCut), opt.depthCut);
    end
end

K = K(K.keep, :);

% ── optional time series ─────────────────────────────────────────────────
if opt.series
    z = cell(height(K),1); v = cell(height(K),1); t = cell(height(K),1);
    for i = 1:height(K)
        [z{i}, v{i}, t{i}] = local_series(K.kinPath(i));
    end
    K.z = z; K.v = v; K.t_s = t;
end

% ── report ───────────────────────────────────────────────────────────────
if opt.verbose
    fprintf('\n--- loaded set ---\n');
    DR = K(~K.isZeroDrop,:);  D0 = K(K.isZeroDrop,:);
    fprintf('  drop-height trials : %3d over %d heights\n', ...
        height(DR), numel(unique(DR.dropHeight_mm)));
    fprintf('  zero-drop trials   : %3d\n', height(D0));
    cs = unique(K.condition);
    for c = 1:numel(cs)
        fprintf('    %-16s %3d\n', cs(c), sum(K.condition==cs(c)));
    end
    if ~isempty(D0)
        fprintf('\n  MEASURED d0 (h = 0, released from contact):\n');
        fprintf('    %.3f +/- %.3f cm   (n = %d, range %.3f - %.3f)\n', ...
            mean(D0.d_final_cm), std(D0.d_final_cm), height(D0), ...
            min(D0.d_final_cm), max(D0.d_final_cm));
        fprintf(['    Ambroso et al. (2005) d+ measured directly rather than\n' ...
                 '    inverted from the drop trials -- compare against the fitted\n' ...
                 '    d0 that depth_scaling(...,''form'',''ambroso'') reports.\n']);
    end
    fprintf('\n');
end
end

% ─────────────────────────────────────────────────────────────────────────
function m = local_model_from_tag(tag)
% Tight and Wide always carry an explicit suffix; everything else is Default.
if endsWith(tag,"_tight"),   m = "Tight";
elseif endsWith(tag,"_wide"), m = "Wide";
else,                         m = "Default";
end
end

function tf = local_isglitch(kinPath)
%LOCAL_ISGLITCH  |v| should fall monotonically from v0 to zero between impact
%   and stop. An increase above 15% of v0 marks a tracking artefact.
tf = false;
if ~isfile(kinPath), return; end
try
    S = load(kinPath,'kin');
catch
    return
end
k = S.kin;
if ~all(isfield(k,{'v','impact_index','stopFrame','v0_cm_s'})), return; end
i1 = max(1, k.impact_index);
i2 = min(numel(k.v), k.stopFrame);
if i2 <= i1 + 1, return; end
seg = abs(k.v(i1:i2));
seg = seg(isfinite(seg));
if numel(seg) < 3, return; end
tf = any(diff(seg) > 0.15*abs(k.v0_cm_s));
end

function [z, v, t] = local_series(kinPath)
z = []; v = []; t = [];
if ~isfile(kinPath), return; end
try
    S = load(kinPath,'kin');
catch
    return
end
k = S.kin;
if isfield(k,'z_cm'),    z = k.z_cm;    elseif isfield(k,'z'), z = k.z; end
if isfield(k,'v_cm_s'),  v = k.v_cm_s;  elseif isfield(k,'v'), v = k.v; end
if isfield(k,'t_s'),     t = k.t_s;     end
end

function v = local_num(T,name)
v = NaN;
if ismember(name,T.Properties.VariableNames)
    x = T.(name)(1);
    if isnumeric(x)||islogical(x), v = double(x); else, v = str2double(string(x)); end
end
end

function s = local_str(T,name)
s = "";
if ismember(name,T.Properties.VariableNames)
    s = strtrim(string(T.(name)(1)));
end
end
