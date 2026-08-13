function OUT = clear_model_kinematics(varargin)
%CLEAR_MODEL_KINEMATICS  Remove kinematics outputs under the model folders.
%
%   DRY RUN BY DEFAULT. Nothing is deleted unless 'Confirm' is true, and even
%   then only the files listed by the dry run.
%
%   SCOPE. Matched on the TRIAL TAG, not on the folder, because the validated
%   269-trial Default dataset lives in the same Default Model tree as the new
%   d0 trials. Only these are matched:
%       *_full_tight     all Tight Model trials
%       *_full_wide      all Wide Model trials
%       0mm_*_full_default   Default Model, ZERO-DROP ONLY
%   The 269 validated Default trials carry no _default suffix, so they cannot
%   match any of these patterns and are structurally out of reach. CHIN and
%   GB/shallow are outside the base path entirely.
%
%   WHAT IT REMOVES. Only kinematics outputs: *_kin.mat, *_kin_scalars.csv and
%   the per-trial kinematics logs. Tracks, detections and exported frames are
%   left alone, so Stage B can regenerate everything in minutes.
%
%   USAGE
%       OUT = clear_model_kinematics;                       % dry run
%       OUT = clear_model_kinematics('Confirm', true);      % delete
%       OUT = clear_model_kinematics('Models', "Tight Model");
%       OUT = clear_model_kinematics('Since', datetime('today'));
%
%   OPTIONS
%       'Root'     default D:\ME_GRANULAB\JerboaImpact
%       'Models'   default ["Default Model","Tight Model","Wide Model"]
%       'Since'    only files modified on or after this datetime ([] = all)
%       'Confirm'  actually delete (default false)

opt.Root    = 'D:\ME_GRANULAB\JerboaImpact';
opt.Models  = ["Default Model","Tight Model","Wide Model"];
opt.Since   = [];
opt.Confirm = false;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

base = fullfile(opt.Root,'03_RESULTS','GB','Batch 5');
fprintf('\n=== clear_model_kinematics ===\n');
fprintf('base   : %s\n', base);
fprintf('models : %s\n', strjoin(cellstr(opt.Models), ', '));
if ~isempty(opt.Since)
    fprintf('since  : %s\n', string(opt.Since));
end
fprintf('mode   : %s\n\n', local_tern(opt.Confirm,'DELETE','DRY RUN'));

rows = {};
% Tag patterns, applied to the FILENAME. This is what protects the 269: they
% are named e.g. 165mm_T03_full_kin.mat, with no model suffix, so no pattern
% below can reach them.
pats = strings(0);
if any(strcmpi(opt.Models,"Tight Model")),   pats(end+1) = "_full_tight";   end
if any(strcmpi(opt.Models,"Wide Model")),    pats(end+1) = "_full_wide";    end
zeroDefault = any(strcmpi(opt.Models,"Default Model"));

F = [dir(fullfile(base,'**','*_kin.mat'));
     dir(fullfile(base,'**','*_kin_scalars.csv'));
     dir(fullfile(base,'**','kinematics','*.txt'))];
F = F(~[F.isdir]);

for i = 1:numel(F)
    nm = string(F(i).name);
    hit = ""; 
    if any(arrayfun(@(p) contains(nm,p), pats))
        if contains(nm,"_full_tight"), hit = "Tight Model"; else, hit = "Wide Model"; end
    elseif zeroDefault && startsWith(nm,"0mm_") && contains(nm,"_full_default")
        hit = "Default Model (d0)";
    end
    if hit == "", continue; end
    dn = datetime(F(i).datenum, 'ConvertFrom','datenum');
    if ~isempty(opt.Since) && dn < opt.Since, continue; end
    rows{end+1} = table(hit, nm, string(F(i).folder), dn, F(i).bytes, ...
        'VariableNames',{'model','file','folder','modified','bytes'}); %#ok<AGROW>
end

if isempty(rows)
    fprintf('Nothing matches. Nothing to do.\n\n');
    OUT = table(); return
end
OUT = vertcat(rows{:});

% ── what would go ────────────────────────────────────────────────────────
fprintf('--- files matched ---\n');
for mm = unique(OUT.model,'stable')'
    s = OUT(OUT.model == mm, :);
    nk = sum(endsWith(s.file,'_kin.mat'));
    nc = sum(endsWith(s.file,'_kin_scalars.csv'));
    no = height(s) - nk - nc;
    fprintf('  %-14s %4d files  (%d kin.mat, %d scalars.csv, %d logs)  %.1f MB\n', ...
        mm, height(s), nk, nc, no, sum(s.bytes)/1e6);
    fprintf('  %-14s modified %s  to  %s\n', '', ...
        string(min(s.modified)), string(max(s.modified)));
end
fprintf('  %-14s %4d files total, %.1f MB\n', 'TOTAL', height(OUT), sum(OUT.bytes)/1e6);

% ── what is deliberately NOT touched ─────────────────────────────────────
fprintf('\n--- untouched (for reassurance) ---\n');
allK = dir(fullfile(base,'**','*_kin.mat'));
allK = allK(~[allK.isdir]);
nm   = string({allK.name}');
spared = ~(contains(nm,"_full_tight") | contains(nm,"_full_wide") | ...
           (startsWith(nm,"0mm_") & contains(nm,"_full_default")));
fprintf('  _kin.mat kept under Batch 5         : %d\n', sum(spared));
if any(spared)
    ex = nm(spared);
    fprintf('    e.g. %s\n', ex(1:min(3,numel(ex))));
end
nonZeroDef = contains(nm,"_full_default") & ~startsWith(nm,"0mm_");
fprintf('  non-zero-drop *_default kin.mat      : %d  (kept)\n', sum(nonZeroDef));
allT = dir(fullfile(base,'**','*_tracks.mat'));
fprintf('  _tracks.mat under Batch 5           : %d  (kept -- Stage B regenerates\n', numel(allT));
fprintf('                                          kinematics from these)\n');

% ── act ──────────────────────────────────────────────────────────────────
if ~opt.Confirm
    fprintf(['\nDRY RUN. Nothing deleted. Re-run with ''Confirm'', true to remove\n' ...
             'exactly the %d files listed above.\n\n'], height(OUT));
    return
end

fprintf('\nDeleting %d files...\n', height(OUT));
nDel = 0; nErr = 0;
for i = 1:height(OUT)
    p = fullfile(OUT.folder(i), OUT.file(i));
    try
        delete(p); nDel = nDel + 1;
    catch ME
        nErr = nErr + 1;
        fprintf('  could not delete %s (%s)\n', OUT.file(i), ME.message);
    end
end
fprintf('deleted %d, failed %d\n', nDel, nErr);

% remove now-empty kinematics folders, but nothing else
K = dir(fullfile(base,'**','kinematics'));
K = K([K.isdir]);
nRm = 0;
for i = 1:numel(K)
    p = fullfile(K(i).folder, K(i).name);
    e = dir(p); e = e(~ismember({e.name},{'.','..'}));
    if isempty(e), rmdir(p); nRm = nRm + 1; end
end
fprintf('removed %d empty kinematics folders\n\n', nRm);
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
