function [DROP, D0, ALL] = load_default_gb(varargin)
%LOAD_DEFAULT_GB  Cleaned Default-geometry GB/full dataset. READ-ONLY.
%
%   Single source of truth for which Default GB/full trials enter the
%   depth-scaling analysis. Both ambroso_collapse and depth_scaling_literature
%   should be driven from this, so the inclusion rule lives in exactly one
%   place and cannot drift between scripts.
%
%   [DROP, D0, ALL] = LOAD_DEFAULT_GB
%       DROP  drop-height trials (h > 0) that survived every rule
%       D0    zero-drop trials (h = 0) that survived the manual list
%       ALL   every candidate considered, with a keep flag and the reason
%             it was dropped -- so an exclusion can always be traced
%
%   SELECTION
%     material/container : GB/full only
%     geometry          : Default only. Two tag conventions coexist --
%                         the validated drop trials carry no model suffix
%                         (165mm_T03_full) while the newer zero-drop trials
%                         do (0mm_T01_full_default). Tight and Wide are
%                         excluded by their _tight / _wide suffixes.
%
%   EXCLUSION RULES
%     1. The manual list below, applied at every height.
%     2. For h > 0 ONLY: NaNDEPTH (non-finite d_final_cm) and GLITCH.
%     3. For h = 0: manual only. No automatic flags, because the criteria
%        below assume a fall and a deceleration that a released-from-rest
%        trial does not have.
%
%   GLITCH is recomputed here from kin.v rather than read from an audit
%   table: |v| should fall monotonically from v0 to zero between impact and
%   stop, so an increase greater than 15% of v0 marks a tracking artefact.
%   The per-rule counts are printed so the outcome can be checked.
%
%   NOTE. The 4.00 cm model-height cutoff that the earlier pipeline applied
%   automatically is NOT applied here, because it was not in the rule set
%   given for this dataset. Pass 'DepthCut' to reinstate it.
%
%   OPTIONS
%       'Root'      default D:\ME_GRANULAB\JerboaImpact
%       'DepthCut'  exclude h > 0 trials deeper than this, cm ([] = off)
%       'Verbose'   print the per-rule breakdown (default true)

opt.Root     = 'D:\ME_GRANULAB\JerboaImpact';
opt.DepthCut = [];
opt.Verbose  = true;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

% ── manual exclusions: the single place this list is written ─────────────
manualExclude = [
    "0mm_T03_full_default"
    "0mm_T05_full_default"
    "0mm_T06_full_default"
    "0mm_T08_full_default"
    "0mm_T09_full_default"
    "25mm_T04_full"
    "45mm_T03_full"
    "65mm_T02_full"
    "85mm_T10_full"
    "205mm_T06_full"
    "205mm_T07_full"
    "205mm_T10_full"
    "225mm_T01_full"
    "245mm_T08_full"
    "245mm_T10_full"
    "285mm_T02_full"
    "285mm_T08_full"
    "305mm_T06_full"
    "325mm_T01_full"
    "325mm_T03_full"
    "325mm_T07_full"
    "325mm_T10_full"
    "345mm_T02_full"
    "345mm_T03_full"
    "345mm_T04_full"
    "345mm_T05_full"
    "345mm_T08_full"
    "345mm_T10_full"
    "365mm_T02_full"
    "365mm_T05_full"
    "365mm_T07_full"
    "365mm_T10_full" ];

fprintf('\n=== load_default_gb ===\n');

% ── gather candidates ────────────────────────────────────────────────────
F = dir(fullfile(opt.Root,'03_RESULTS','GB','**','*_kin_scalars.csv'));
F = F(~[F.isdir]);
if isempty(F)
    error('load_default_gb:noFiles','No GB *_kin_scalars.csv under %s', opt.Root);
end

n = numel(F);
tag=strings(n,1); cond=strings(n,1); h=nan(n,1); d=nan(n,1); v0=nan(n,1);
ts=nan(n,1); phi=nan(n,1); kp=strings(n,1);
for i = 1:n
    T = readtable(fullfile(F(i).folder,F(i).name));
    if height(T) < 1, continue; end
    tag(i)  = string(erase(F(i).name,'_kin_scalars.csv'));
    cond(i) = local_str(T,'condition');
    d(i)    = local_num(T,'d_final_cm');
    v0(i)   = abs(local_num(T,'v0_cm_s'));
    ts(i)   = local_num(T,'t_stop_s');
    phi(i)  = local_num(T,'phi');
    hh      = local_num(T,'dropHeight_true_mm');
    if ~isfinite(hh), hh = local_num(T,'dropHeight_mm'); end
    h(i)    = hh;
    kp(i)   = string(strrep(fullfile(F(i).folder,F(i).name), ...
                            '_kin_scalars.csv','_kin.mat'));
end

% ── geometry + condition filter ──────────────────────────────────────────
isGBfull  = cond == "GB/full";
isTight   = endsWith(tag,"_tight");
isWide    = endsWith(tag,"_wide");
isDefault = (endsWith(tag,"_full") | endsWith(tag,"_full_default")) & ~isTight & ~isWide;
cand      = tag ~= "" & isGBfull & isDefault;

fprintf('GB kin_scalars found      : %d\n', sum(tag ~= ""));
fprintf('  GB/full & Default       : %d\n', sum(cand));
fprintf('  excluded, other geometry: %d tight, %d wide\n', ...
        sum(isTight & isGBfull), sum(isWide & isGBfull));
fprintf('  excluded, other condition: %d\n', sum(tag ~= "" & ~isGBfull));

ALL = table(tag(cand), cond(cand), h(cand), d(cand), v0(cand), ts(cand), ...
            phi(cand), kp(cand), ...
    'VariableNames',{'trialTag','condition','dropHeight_mm','d_final_cm', ...
                     'v0_cm_s','t_stop_s','phi','kinPath'});
ALL.isZeroDrop = ALL.dropHeight_mm == 0;

% ── rules ────────────────────────────────────────────────────────────────
ALL.exManual   = ismember(ALL.trialTag, manualExclude);
ALL.exNaNDepth = ~ALL.isZeroDrop & ~isfinite(ALL.d_final_cm);

ALL.exGlitch = false(height(ALL),1);
for i = 1:height(ALL)
    if ALL.isZeroDrop(i) || ALL.exManual(i) || ALL.exNaNDepth(i), continue; end
    ALL.exGlitch(i) = local_isglitch(ALL.kinPath(i));
end

ALL.exDepthCut = false(height(ALL),1);
if ~isempty(opt.DepthCut)
    ALL.exDepthCut = ~ALL.isZeroDrop & ALL.d_final_cm > opt.DepthCut;
end

ALL.keep = ~(ALL.exManual | ALL.exNaNDepth | ALL.exGlitch | ALL.exDepthCut);
ALL.reason = strings(height(ALL),1);
ALL.reason(ALL.exManual)   = "manual";
ALL.reason(ALL.exNaNDepth) = "NaNDEPTH";
ALL.reason(ALL.exGlitch)   = "GLITCH";
ALL.reason(ALL.exDepthCut) = "depth cut";

DROP = ALL(ALL.keep & ~ALL.isZeroDrop, :);
D0   = ALL(ALL.keep &  ALL.isZeroDrop, :);

% ── report ───────────────────────────────────────────────────────────────
if opt.Verbose
    fprintf('\n--- exclusions applied ---\n');
    fprintf('  manual   : %2d of %d listed  (%d listed tags matched no file)\n', ...
        sum(ALL.exManual), numel(manualExclude), ...
        sum(~ismember(manualExclude, ALL.trialTag)));
    miss = manualExclude(~ismember(manualExclude, ALL.trialTag));
    if ~isempty(miss), fprintf('             %s\n', miss); end
    fprintf('  NaNDEPTH : %2d  (h > 0 only)\n', sum(ALL.exNaNDepth));
    fprintf('  GLITCH   : %2d  (h > 0 only)\n', sum(ALL.exGlitch));
    if ~isempty(opt.DepthCut)
        fprintf('  depth cut: %2d  (> %.2f cm)\n', sum(ALL.exDepthCut), opt.DepthCut);
    end

    fprintf('\n--- cleaned dataset ---\n');
    fprintf('  drop-height trials : %3d over %d heights\n', ...
        height(DROP), numel(unique(DROP.dropHeight_mm)));
    fprintf('  zero-drop trials   : %3d\n', height(D0));

    fprintf('\n  replicates per height:\n');
    hs = unique(DROP.dropHeight_mm);
    for k = 1:numel(hs)
        fprintf('    %4g mm : %2d\n', hs(k), sum(DROP.dropHeight_mm == hs(k)));
    end

    if ~isempty(D0)
        fprintf('\n  MEASURED d0 (h = 0, released from contact):\n');
        fprintf('    %.3f +/- %.3f cm   (n = %d, range %.3f - %.3f)\n', ...
            mean(D0.d_final_cm), std(D0.d_final_cm), height(D0), ...
            min(D0.d_final_cm), max(D0.d_final_cm));
        fprintf(['    This is Ambroso et al. (2005) d+ measured directly rather\n' ...
                 '    than inverted from the drop trials -- compare it against the\n' ...
                 '    fitted d0 that ambroso_collapse reports.\n']);
    end
    fprintf('\n');
end
end

% ─────────────────────────────────────────────────────────────────────────
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

function v = local_num(T,name)
v = NaN;
if ismember(name,T.Properties.VariableNames)
    x = T.(name)(1);
    if isnumeric(x)||islogical(x), v = double(x); else, v = str2double(string(x)); end
end
end
function s = local_str(T,name)
s = ""; if ismember(name,T.Properties.VariableNames), s = strtrim(string(T.(name)(1))); end
end
