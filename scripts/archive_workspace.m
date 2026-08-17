function R = archive_workspace(varargin)
%ARCHIVE_WORKSPACE  Move processed model subtrees out to the archive and zip them.
%
%   MAINTAINER TOOL. DRY RUN BY DEFAULT -- prints the whole plan and touches
%   nothing. Pass 'Confirm',true to actually move, zip, verify and delete.
%
%       archive_workspace                                  % plan only
%       archive_workspace('Root','D:\ME_GRANULAB\JerboaImpact')
%       archive_workspace('Confirm',true)                  % do it
%
%   WHAT IT TOUCHES
%     Only these three trees under Root:
%       01_FRAMES  02_SAVED_DETECTIONS  03_RESULTS
%     and within each, only GB\<batch>\<model> subtrees where <model> is one of
%     the recognised foot models (Default / Tight / Wide Model). Every batch
%     level present is picked up.
%
%   WHAT IT NEVER TOUCHES
%     - raw videos: they live outside these three trees and are never scanned
%     - anything already inside ArchiveRoot: existing archive contents are read
%       for collision checks only, never moved, overwritten or deleted
%     - CHIN, or any GB subtree with no model level (the pre-model layout)
%
%   ORDER OF OPERATIONS (each step completes for all sources before the next)
%     1. Copy every DEFAULT *_kin_scalars.csv flat into
%           <ArchiveRoot>\reference_default_scalars\
%        UNCOMPRESSED and left there permanently. This is the baseline for the
%        old-vs-new impactDistPx (-360) comparison, so it must stay readable
%        without unzipping anything. Filenames are already unique per trial;
%        a collision is reported and the file skipped rather than overwritten.
%     2. Move each source subtree, hierarchy preserved, to
%           <ArchiveRoot>\<MODEL>\<BATCH>\<DATA_TYPE>\...
%        where DATA_TYPE is frames | detections | results.
%     3. Zip one archive per <MODEL>\<BATCH>, covering all three data types.
%     4. VERIFY by file count -- entries in the zip against the files that were
%        moved -- and only then delete the uncompressed copy. On any mismatch
%        BOTH are kept and a warning is issued: a zip that cannot be shown to
%        be complete is not allowed to become the only copy.
%     5. Write <ArchiveRoot>\ARCHIVE_MANIFEST.txt.
%
%   OPTIONS (name-value)
%       'Root'        workspace root       (default D:\ME_GRANULAB\JerboaImpact)
%       'FramesRoot'  where 01_FRAMES lives; frames may sit on a local scratch
%                     disk (default: $JERBOA_FRAMES_ROOT, else Root -- the same
%                     rule process_trial uses)
%       'ArchiveRoot' destination         (default D:\ME_GRANULAB\_ARCHIVE)
%       'Models'      recognised model folder names
%                     (default ["Default Model","Tight Model","Wide Model"])
%       'Confirm'     false (default) = dry run. true = act.
%       'Zip'         true (default). false moves but leaves it uncompressed,
%                     which also means nothing is deleted.
%
%   Returns R, the per-subtree plan/result table, so a dry run can be inspected
%   in the workspace before committing.

opt.Root        = 'D:\ME_GRANULAB\JerboaImpact';
opt.FramesRoot  = '';
opt.ArchiveRoot = 'D:\ME_GRANULAB\_ARCHIVE';
opt.Models      = ["Default Model","Tight Model","Wide Model"];
opt.Confirm     = false;
opt.Zip         = true;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

if isempty(opt.FramesRoot)
    opt.FramesRoot = getenv('JERBOA_FRAMES_ROOT');
    if isempty(opt.FramesRoot), opt.FramesRoot = opt.Root; end
end
opt.Models = string(opt.Models);

DATA_TYPES = { ...
    '01_FRAMES',           'frames',     opt.FramesRoot ; ...
    '02_SAVED_DETECTIONS', 'detections', opt.Root       ; ...
    '03_RESULTS',          'results',    opt.Root       };

fprintf('\n=== archive_workspace ===\n');
fprintf('  Root        : %s\n', opt.Root);
fprintf('  FramesRoot  : %s\n', opt.FramesRoot);
fprintf('  ArchiveRoot : %s\n', opt.ArchiveRoot);
fprintf('  Models      : %s\n', strjoin(cellstr(opt.Models), ', '));
fprintf('  Zip         : %s\n', local_tern(opt.Zip,'yes','no'));
fprintf('  MODE        : %s\n\n', local_tern(opt.Confirm, ...
        '*** CONFIRM -- files WILL be moved and deleted ***', ...
        'DRY RUN (nothing will be touched)'));

local_guard_paths(opt);

% ── 1) discover source subtrees ──────────────────────────────────────────
rows = {};
for d = 1:size(DATA_TYPES,1)
    typeDir  = DATA_TYPES{d,1};
    typeName = DATA_TYPES{d,2};
    typeRoot = fullfile(DATA_TYPES{d,3}, typeDir);
    if ~isfolder(typeRoot)
        fprintf('  (absent, skipped) %s\n', typeRoot);
        continue
    end
    gbRoot = fullfile(typeRoot, 'GB');
    if ~isfolder(gbRoot)
        fprintf('  (no GB tree)      %s\n', typeRoot);
        continue
    end
    batches = local_subdirs(gbRoot);
    for b = 1:numel(batches)
        batchDir = batches{b};
        models   = local_subdirs(fullfile(gbRoot, batchDir));
        for mi = 1:numel(models)
            if ~any(strcmpi(models{mi}, opt.Models)), continue; end
            src = fullfile(gbRoot, batchDir, models{mi});
            [nFiles, nBytes] = local_tree_size(src);
            rows{end+1} = table(string(models{mi}), string(batchDir), ...
                string(typeName), string(src), nFiles, nBytes/2^30, ...
                'VariableNames', {'model','batch','dataType','source', ...
                                  'nFiles','GB'}); %#ok<AGROW>
        end
    end
end

if isempty(rows)
    fprintf(['\nNothing to archive: no GB\\<batch>\\<model> subtree found under\n' ...
             '01_FRAMES, 02_SAVED_DETECTIONS or 03_RESULTS.\n' ...
             'Check Root/FramesRoot, and that the model folders are named %s.\n\n'], ...
             strjoin(cellstr(opt.Models), ' / '));
    R = table();
    return
end
R = vertcat(rows{:});
R.dest = strings(height(R),1);
for i = 1:height(R)
    R.dest(i) = string(fullfile(opt.ArchiveRoot, char(R.model(i)), ...
                                char(R.batch(i)), char(R.dataType(i))));
end

% ── 2) plan ──────────────────────────────────────────────────────────────
fprintf('\n--- PLAN ---\n');
pairs = unique(R(:,{'model','batch'}), 'rows');
for p = 1:height(pairs)
    sel = R.model==pairs.model(p) & R.batch==pairs.batch(p);
    fprintf('\n  %s \\ %s   (%d files, %.2f GB)\n', ...
            pairs.model(p), pairs.batch(p), sum(R.nFiles(sel)), sum(R.GB(sel)));
    idx = find(sel).';
    for i = idx
        fprintf('    %-11s %6d files  %7.2f GB\n', R.dataType(i), R.nFiles(i), R.GB(i));
        fprintf('        from %s\n', R.source(i));
        fprintf('        to   %s\n', R.dest(i));
    end
    if opt.Zip
        fprintf('    zip -> %s\n', local_zip_path(opt, pairs.model(p), pairs.batch(p)));
    end
end
fprintf('\n  TOTAL: %d files, %.2f GB across %d subtree(s)\n', ...
        sum(R.nFiles), sum(R.GB), height(R));

refCsvs = local_find_default_scalars(R, opt);
refDir  = fullfile(opt.ArchiveRoot, 'reference_default_scalars');
fprintf('\n  reference copy (uncompressed, kept permanently):\n');
fprintf('    %d Default *_kin_scalars.csv -> %s\n', numel(refCsvs), refDir);

if ~opt.Confirm
    fprintf(['\nDRY RUN -- nothing was touched.\n' ...
             'Re-run with ''Confirm'',true to move, zip, verify and delete.\n\n']);
    return
end

% ── 3) reference scalars copy (BEFORE anything moves) ────────────────────
fprintf('\n--- COPYING REFERENCE SCALARS ---\n');
if ~isfolder(refDir), mkdir(refDir); end
nRef = 0; nRefSkip = 0;
for i = 1:numel(refCsvs)
    [~, nm, ex] = fileparts(refCsvs{i});
    target = fullfile(refDir, [nm ex]);
    if isfile(target)
        warning('archive_workspace:refCollision', ...
            'Reference copy already exists, not overwritten: %s', target);
        nRefSkip = nRefSkip + 1;
        continue
    end
    copyfile(refCsvs{i}, target);       % COPY, not move: the original still
    nRef = nRef + 1;                    % travels with its subtree
end
fprintf('  copied %d, skipped %d (already present)\n', nRef, nRefSkip);

% ── 4) move ──────────────────────────────────────────────────────────────
fprintf('\n--- MOVING ---\n');
R.moved      = false(height(R),1);
R.movedFiles = zeros(height(R),1);
for i = 1:height(R)
    src = char(R.source(i));  dst = char(R.dest(i));
    if isfolder(dst)
        warning('archive_workspace:destExists', ...
            ['Destination already exists, subtree NOT moved: %s\n' ...
             'Existing archive contents are never overwritten. Resolve by hand.'], dst);
        continue
    end
    parent = fileparts(dst);
    if ~isfolder(parent), mkdir(parent); end
    [ok, msg] = movefile(src, dst);
    if ~ok
        warning('archive_workspace:moveFailed','%s -> %s : %s', src, dst, msg);
        continue
    end
    R.moved(i)      = true;
    R.movedFiles(i) = local_tree_size(dst);
    fprintf('  moved %-11s %6d files  %s\n', R.dataType(i), R.movedFiles(i), dst);
end

% ── 5) zip, verify, then delete ──────────────────────────────────────────
R.zipPath = strings(height(R),1);
R.verified = false(height(R),1);
zipRows = {};
if opt.Zip
    fprintf('\n--- ZIPPING AND VERIFYING ---\n');
    for p = 1:height(pairs)
        sel = find(R.model==pairs.model(p) & R.batch==pairs.batch(p) & R.moved).';
        if isempty(sel), continue; end
        pairDir = fullfile(opt.ArchiveRoot, char(pairs.model(p)), char(pairs.batch(p)));
        zipPath = local_zip_path(opt, pairs.model(p), pairs.batch(p));
        expected = sum(R.movedFiles(sel));

        if isfile(zipPath)
            warning('archive_workspace:zipExists', ...
                'Zip already exists, not overwritten: %s', zipPath);
            continue
        end

        fprintf('  %s \\ %s -> %s\n', pairs.model(p), pairs.batch(p), zipPath);
        try
            % The zip is written at pairDir\<name>.zip while the entries are the
            % frames/detections/results subfolders, so the archive is never
            % inside its own input and cannot capture itself mid-write.
            zip(zipPath, cellstr(R.dataType(sel)), pairDir);
        catch ME
            warning('archive_workspace:zipFailed', ...
                'zip failed for %s (%s). Uncompressed copy KEPT.', pairDir, ME.message);
            continue
        end

        nInZip = local_zip_entry_count(zipPath);
        okVerify = isfinite(nInZip) && nInZip == expected;
        fprintf('    verify: %d files moved, %d entries in zip -> %s\n', ...
                expected, nInZip, local_tern(okVerify,'MATCH','MISMATCH'));

        if okVerify
            R.verified(sel) = true;
            R.zipPath(sel)  = string(zipPath);
            for i = sel
                rmdir(char(R.dest(i)), 's');       % only now is it safe
            end
            fprintf('    uncompressed copy deleted\n');
        else
            R.zipPath(sel) = string(zipPath);
            warning('archive_workspace:verifyFailed', ...
                ['File-count mismatch for %s (moved %d, zip %d). BOTH the zip ' ...
                 'and the uncompressed copy have been KEPT. Do not delete either ' ...
                 'until the difference is understood.'], pairDir, expected, nInZip);
        end
        zipRows{end+1} = struct('model',pairs.model(p),'batch',pairs.batch(p), ...
            'zipPath',string(zipPath),'expected',expected,'inZip',nInZip, ...
            'verified',okVerify); %#ok<AGROW>
    end
end

% ── 6) manifest ──────────────────────────────────────────────────────────
local_write_manifest(opt, R, zipRows, nRef, nRefSkip, refDir);

fprintf('\n=== ARCHIVE COMPLETE ===\n');
fprintf('  subtrees moved   : %d of %d\n', sum(R.moved), height(R));
if opt.Zip
    fprintf('  zips verified    : %d\n', numel(unique(R.zipPath(R.verified))));
    nUnver = sum(R.moved & ~R.verified);
    if nUnver > 0
        fprintf('  UNVERIFIED       : %d subtree(s) still uncompressed on disk\n', nUnver);
    end
end
fprintf('  manifest         : %s\n\n', fullfile(opt.ArchiveRoot,'ARCHIVE_MANIFEST.txt'));
end

% ═════════════════════════════════════════════════════════════════════════
function local_guard_paths(opt)
%LOCAL_GUARD_PATHS  Refuse configurations that could archive into the workspace.
    a = local_norm(opt.ArchiveRoot);
    for r = {local_norm(opt.Root), local_norm(opt.FramesRoot)}
        base = r{1};
        for t = {'01_FRAMES','02_SAVED_DETECTIONS','03_RESULTS'}
            live = local_norm(fullfile(base, t{1}));
            if startsWith(a, live)
                error('archive_workspace:archiveInsideWorkspace', ...
                    ['ArchiveRoot (%s) is inside the live tree %s. The archive ' ...
                     'must sit outside the folders being archived, or the move ' ...
                     'would consume its own destination.'], opt.ArchiveRoot, live);
            end
        end
    end
end

function p = local_norm(p)
    p = lower(strrep(char(p), '/', filesep));
    if endsWith(p, filesep), p = p(1:end-1); end
end

function d = local_subdirs(root)
    d = {};
    if ~isfolder(root), return; end
    L = dir(root);
    L = L([L.isdir]);
    d = {L(~ismember({L.name},{'.','..'})).name};
end

function [nFiles, nBytes] = local_tree_size(root)
%LOCAL_TREE_SIZE  Recursive file count and total bytes. Files only, no dirs.
    nFiles = 0; nBytes = 0;
    if ~isfolder(root), return; end
    L = dir(fullfile(root, '**', '*'));
    L = L(~[L.isdir]);
    nFiles = numel(L);
    nBytes = sum([L.bytes]);
end

function csvs = local_find_default_scalars(R, opt) %#ok<INUSD>
%LOCAL_FIND_DEFAULT_SCALARS  Every *_kin_scalars.csv under a Default results
%   subtree that is about to be archived. Results only -- scalars live there.
    csvs = {};
    sel = find(startsWith(lower(R.model), 'default') & R.dataType == "results").';
    for i = sel
        L = dir(fullfile(char(R.source(i)), '**', '*_kin_scalars.csv'));
        L = L(~[L.isdir]);
        for k = 1:numel(L)
            csvs{end+1} = fullfile(L(k).folder, L(k).name); %#ok<AGROW>
        end
    end
end

function p = local_zip_path(opt, model, batch)
    safeModel = matlab.lang.makeValidName(char(model));
    safeBatch = matlab.lang.makeValidName(char(batch));
    p = fullfile(opt.ArchiveRoot, char(model), char(batch), ...
                 sprintf('%s_%s.zip', safeModel, safeBatch));
end

function n = local_zip_entry_count(zipPath)
%LOCAL_ZIP_ENTRY_COUNT  Count FILE entries by reading the written zip back.
%   Reads the artefact itself rather than trusting the zip() call, since the
%   whole point of the check is to prove the zip is complete before anything is
%   deleted. Directory entries are excluded so the count is comparable with a
%   files-only tree count. Returns NaN if the zip cannot be read, which the
%   caller treats as a failed verification.
    n = NaN;
    try
        zf = java.util.zip.ZipFile(java.io.File(zipPath));
        c  = 0;
        en = zf.entries();
        while en.hasMoreElements()
            e = en.nextElement();
            if ~e.isDirectory(), c = c + 1; end
        end
        zf.close();
        n = c;
    catch
        % leave NaN
    end
end

function local_write_manifest(opt, R, zipRows, nRef, nRefSkip, refDir)
    mp = fullfile(opt.ArchiveRoot, 'ARCHIVE_MANIFEST.txt');
    if ~isfolder(opt.ArchiveRoot), mkdir(opt.ArchiveRoot); end
    fid = fopen(mp, 'a');            % APPEND: earlier archive runs are history
    if fid < 0
        warning('archive_workspace:manifestFailed','Could not write %s', mp);
        return
    end
    fprintf(fid, '\n======================================================================\n');
    fprintf(fid, ' ARCHIVE RUN  %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '======================================================================\n');
    fprintf(fid, 'Root        : %s\n', opt.Root);
    fprintf(fid, 'FramesRoot  : %s\n', opt.FramesRoot);
    fprintf(fid, 'ArchiveRoot : %s\n\n', opt.ArchiveRoot);

    fprintf(fid, 'REFERENCE SCALARS (uncompressed, kept permanently)\n');
    fprintf(fid, '  %s\n', refDir);
    fprintf(fid, '  copied %d, skipped %d (already present)\n\n', nRef, nRefSkip);

    fprintf(fid, 'SUBTREES\n');
    for i = 1:height(R)
        fprintf(fid, '  [%s] %s \\ %s \\ %s\n', ...
            local_tern(R.moved(i),'MOVED','NOT MOVED'), ...
            R.model(i), R.batch(i), R.dataType(i));
        fprintf(fid, '      source : %s\n', R.source(i));
        fprintf(fid, '      dest   : %s\n', R.dest(i));
        fprintf(fid, '      files  : %d source, %d moved   (%.2f GB)\n', ...
            R.nFiles(i), R.movedFiles(i), R.GB(i));
    end

    fprintf(fid, '\nZIPS AND VERIFICATION\n');
    if isempty(zipRows)
        fprintf(fid, '  (none)\n');
    else
        for z = 1:numel(zipRows)
            zr = zipRows{z};
            fprintf(fid, '  %s \\ %s\n', zr.model, zr.batch);
            fprintf(fid, '      zip      : %s\n', zr.zipPath);
            fprintf(fid, '      expected : %d files\n', zr.expected);
            fprintf(fid, '      in zip   : %d entries\n', zr.inZip);
            if zr.verified
                fprintf(fid, '      VERIFIED : uncompressed copy deleted\n');
            else
                fprintf(fid, '      MISMATCH : BOTH copies kept, nothing deleted\n');
            end
        end
    end

    fprintf(fid, '\nTOTAL: %d files, %.2f GB across %d subtree(s)\n', ...
        sum(R.nFiles), sum(R.GB), height(R));
    fprintf(fid, '======================================================================\n');
    fclose(fid);
    fprintf('\n  manifest appended: %s\n', mp);
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
