function OUT = scan_new_models(varargin)
%SCAN_NEW_MODELS  Dry run over the new Batch 5 GB/full videos. READ-ONLY.
%
%   Writes nothing, processes nothing. Reports what is there, what is already
%   processed, and what will NOT work if run as-is.
%
%   USAGE
%       OUT = scan_new_models;
%       OUT = scan_new_models('Root','F:\ME_GRANULAB\JerboaImpact');
%
%   OPTIONS
%       'RawRoot'  raw video root (default F:\ME_GRANULAB\Test Batches\Batch 5)
%       'Root'     results root   (default F:\ME_GRANULAB\JerboaImpact)
%       'Models'   model folders to scan (default the three below)

opt.RawRoot = 'F:\ME_GRANULAB\Test Batches\Batch 5';
opt.Root    = 'F:\ME_GRANULAB\JerboaImpact';
opt.Models  = ["Default Model","Tight Model","Wide Model"];
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

fprintf('\n=== scan_new_models  (READ-ONLY dry run) ===\n');

% existing outputs, for the "already processed" check
Etr = dir(fullfile(opt.Root,'03_RESULTS','**','*_tracks.mat'));
Ekn = dir(fullfile(opt.Root,'03_RESULTS','**','*_kin.mat'));
trTags  = string(erase({Etr.name},'_tracks.mat'))';
knTags  = string(erase({Ekn.name},'_kin.mat'))';
trPaths = string(fullfile({Etr.folder}',{Etr.name}'));
fprintf('existing outputs: %d tracks, %d kin\n\n', numel(trTags), numel(knTags));

rows = {};
for m = 1:numel(opt.Models)
    model = opt.Models(m);
    fdir  = fullfile(opt.RawRoot, char(model), 'GB', 'full');
    fprintf('--- %s ---\n%s\n', model, fdir);
    if ~isfolder(fdir)
        fprintf('  FOLDER NOT FOUND -- check the path\n\n');
        continue
    end
    V = dir(fullfile(fdir,'*.avi'));
    V = V(~[V.isdir]);
    fprintf('  videos found: %d\n', numel(V));
    if isempty(V), fprintf('\n'); continue; end

    for i = 1:numel(V)
        [~, stem] = fileparts(V(i).name);
        % expected stem: <h>mm_T##
        tok = regexp(stem, '^(\d+)mm_T(\d+)$', 'tokens', 'once');
        if isempty(tok)
            parses = false; hmm = NaN; tnum = NaN;
        else
            parses = true; hmm = str2double(tok{1}); tnum = str2double(tok{2});
        end
        baseTag = string(stem) + "_full";          % what the pipeline builds today
        collide = ismember(baseTag, trTags);       % same tag already on disk

        rows{end+1} = table(model, string(V(i).name), string(stem), parses, ...
            hmm, tnum, baseTag, collide, ismember(baseTag,knTags), ...
            hmm==0, V(i).bytes/1e9, ...
            'VariableNames',{'model','file','stem','parses','dropHeight_mm', ...
                'trialNum','baseTag','tagExistsInResults','hasKin','isZeroDrop', ...
                'sizeGB'}); %#ok<AGROW>
    end
    fprintf('\n');
end
if isempty(rows), error('scan_new_models:none','No videos found.'); end
OUT.files = vertcat(rows{:});

% ------------------------------------------------------------ per model
fprintf('--- per model ---\n');
for m = unique(OUT.files.model,'stable')'
    T = OUT.files(OUT.files.model==m,:);
    fprintf('  %-14s %2d videos | %2d zero-drop | %2d unparseable | %2d tags already in results\n', ...
        m, height(T), sum(T.isZeroDrop), sum(~T.parses), sum(T.tagExistsInResults));
    hs = unique(T.dropHeight_mm(T.parses));
    fprintf('                 heights: %s\n', strjoin(compose('%g',hs'), ', '));
end

% ------------------------------------------------------------ problems
fprintf('\n--- ISSUES THAT MUST BE RESOLVED BEFORE PROCESSING ---\n');

bad = OUT.files(~OUT.files.parses,:);
if ~isempty(bad)
    fprintf('\n1) Filenames that will not parse as <h>mm_T##:\n');
    fprintf('   %-14s %s\n', bad.model', bad.file');
else
    fprintf('\n1) Filename parsing: all OK.\n');
end

col = OUT.files(OUT.files.tagExistsInResults,:);
if ~isempty(col)
    fprintf(['\n2) TAG COLLISIONS -- %d file(s) would produce a trialTag that ALREADY\n' ...
             '   exists in 03_RESULTS. The tag is built from height + trial + container\n' ...
             '   and carries NO model, so Tight/Wide trials collide with Default ones.\n' ...
             '   Processing these as-is would overwrite existing Default Model output.\n'], height(col));
    for i = 1:min(10,height(col))
        j = find(trTags==col.baseTag(i),1);
        fprintf('   %-14s %-16s -> %s\n      already at: %s\n', ...
            col.model(i), col.file(i), col.baseTag(i), trPaths(j));
    end
    if height(col) > 10, fprintf('   ... and %d more\n', height(col)-10); end
else
    fprintf('\n2) Tag collisions: none.\n');
end

z = OUT.files(OUT.files.isZeroDrop,:);
if ~isempty(z)
    fprintf(['\n3) ZERO-DROP TRIALS -- %d file(s) at 0 mm.\n' ...
             '   These are d0 measurements: the foot starts in contact, so there is no\n' ...
             '   free-fall segment and no velocity peak. kd_kinematics locates impact by\n' ...
             '   the impactDistPx trigger refined to the velocity peak, so it will not\n' ...
             '   find a meaningful impact frame here and d_final would come out near zero\n' ...
             '   regardless of the true sink depth. These need a separate measurement\n' ...
             '   (surface reference frame vs. rest frame), not the impact pipeline.\n'], height(z));
    fprintf('   %s\n', z.file');
end

fprintf(['\n4) Model identity is NOT currently carried anywhere: not in the tag, not in\n' ...
         '   meta, and not in the results path (<GB|CHIN>\\Batch 5\\<h>mm_T##\\<container>).\n' ...
         '   A model field and a model level in the output path are needed before any\n' ...
         '   Tight or Wide trial is written.\n']);

fprintf('\nNothing was written. OUT.files returned.\n\n');
end
