function OUT = clear_model_all(varargin)
%CLEAR_MODEL_ALL  Remove ALL processing outputs (Stage A + Stage B) for the
%   Tight/Wide foot models: exported frames, detections, tracks, QA, and
%   kinematics.
%
%   DRY RUN BY DEFAULT. Nothing is deleted unless 'Confirm' is true, and even
%   then only the directories listed by the dry run.
%
%   SCOPE. Matched on the MODEL FOLDER, not on tags, because the model level
%   is a directory in every output tree:
%       <root>/01_FRAMES/GB/Batch 5/<Model>/...            (or under
%           JERBOA_FRAMES_ROOT if that env var was set at processing time --
%           BOTH locations are checked, since some runs may predate it)
%       <root>/02_SAVED_DETECTIONS/GB/Batch 5/<Model>/...
%       <root>/03_RESULTS/GB/Batch 5/<Model>/...
%   Those folders contain ONLY that model's outputs, so removing them erases
%   every trace of its processing and cannot reach anything else.
%
%   'Default Model' is REFUSED unconditionally: the validated 269-trial
%   dataset lives in that tree (see clear_model_kinematics). If Default's
%   zero-drop trials ever need clearing, use clear_model_kinematics for the
%   kinematics and handle tracks by hand, deliberately.
%
%   Top-level 03_RESULTS/_batch_logs entries (video lists, kin logs) are kept
%   as history; they reference nothing after deletion and are harmless.
%
%   USAGE
%       OUT = clear_model_all;                      % dry run, Tight + Wide
%       OUT = clear_model_all('Confirm', true);     % delete
%       OUT = clear_model_all('Models', "Wide Model");
%
%   OPTIONS
%       'Root'        default D:\ME_GRANULAB\JerboaImpact
%       'FramesRoot'  default '' = getenv('JERBOA_FRAMES_ROOT'), else Root
%       'BatchLabel'  default 'Batch 5'
%       'Models'      default ["Tight Model","Wide Model"]
%       'Confirm'     actually delete (default false)

opt.Root       = 'D:\ME_GRANULAB\JerboaImpact';
opt.FramesRoot = '';
opt.BatchLabel = 'Batch 5';
opt.Models     = ["Tight Model","Wide Model"];
opt.Confirm    = false;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

% ── hard guard ────────────────────────────────────────────────────────────
if any(contains(lower(string(opt.Models)), "default"))
    error('clear_model_all:refuseDefault', ...
        ['Default Model is refused: the validated 269-trial dataset lives in ' ...
         'that tree. This script will not touch it.']);
end

if isempty(opt.FramesRoot)
    opt.FramesRoot = getenv('JERBOA_FRAMES_ROOT');
    if isempty(opt.FramesRoot), opt.FramesRoot = opt.Root; end
end

fprintf('\n=== clear_model_all ===\n');
fprintf('root       : %s\n', opt.Root);
fprintf('framesRoot : %s\n', opt.FramesRoot);
fprintf('models     : %s\n', strjoin(cellstr(opt.Models), ', '));
fprintf('mode       : %s\n\n', local_tern(opt.Confirm, 'DELETE', 'DRY RUN'));

% ── enumerate candidate model directories ─────────────────────────────────
% Frames are checked under BOTH FramesRoot and Root (env-var and legacy runs);
% duplicates collapse via unique().
cand = strings(0,1);
for m = 1:numel(opt.Models)
    model = char(opt.Models(m));
    cand(end+1) = string(fullfile(opt.FramesRoot,'01_FRAMES','GB',opt.BatchLabel,model)); %#ok<AGROW>
    cand(end+1) = string(fullfile(opt.Root,      '01_FRAMES','GB',opt.BatchLabel,model)); %#ok<AGROW>
    cand(end+1) = string(fullfile(opt.Root,'02_SAVED_DETECTIONS','GB',opt.BatchLabel,model)); %#ok<AGROW>
    cand(end+1) = string(fullfile(opt.Root,'03_RESULTS',         'GB',opt.BatchLabel,model)); %#ok<AGROW>
end
cand = unique(cand, 'stable');

rows = {};
for i = 1:numel(cand)
    p = char(cand(i));
    if ~isfolder(p)
        fprintf('  absent : %s\n', p);
        continue
    end
    F = dir(fullfile(p,'**','*'));  F = F(~[F.isdir]);
    nB = sum([F.bytes]);
    fprintf('  FOUND  : %s\n           %6d files, %8.1f MB\n', p, numel(F), nB/1e6);
    rows{end+1} = table(string(p), numel(F), nB, 'VariableNames', ...
                        {'dir','nFiles','bytes'}); %#ok<AGROW>
end

if isempty(rows)
    fprintf('\nNothing found to delete.\n\n');
    OUT = table(); return
end
OUT = vertcat(rows{:});
fprintf('\nTOTAL: %d dirs, %d files, %.1f MB\n', ...
        height(OUT), sum(OUT.nFiles), sum(OUT.bytes)/1e6);

% ── delete on confirm ─────────────────────────────────────────────────────
if ~opt.Confirm
    fprintf('\nDRY RUN -- nothing deleted. Re-run with ''Confirm'', true.\n\n');
    return
end
for i = 1:height(OUT)
    p = char(OUT.dir(i));
    [ok, msg] = rmdir(p, 's');
    fprintf('  %s : %s\n', local_tern(ok,'deleted','FAILED'), p);
    if ~ok, warning('clear_model_all:rmdir','%s', msg); end
end
fprintf('\nDone. Stage A can regenerate everything from the raw .avi files.\n\n');
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
