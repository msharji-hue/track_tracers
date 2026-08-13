function write_height_labels(root, varargin)
% WRITE_HEIGHT_LABELS  Drop human-readable warning files into the results tree
% documenting the GB/shallow drop-height label reversal.
%
%   Writes ONLY .txt marker files. Never renames, moves, edits or deletes any
%   data file, folder or trialTag.
%
%   write_height_labels(root)                % write the markers
%   write_height_labels(root,'dryRun',true)  % list what would be written
%   write_height_labels(root,'remove',true)  % delete the markers again
%
%   Creates:
%     03_RESULTS/HEIGHT_LABELS_README.txt          top-level explanation + map
%     .../<any GB shallow container>/_HEIGHT_REVERSED.txt   per-folder marker

    p = inputParser;
    addParameter(p,'dryRun',false,@islogical);
    addParameter(p,'remove',false,@islogical);
    parse(p,varargin{:}); o = p.Results;

    resultsRoot = fullfile(root,'03_RESULTS');
    if ~isfolder(resultsRoot), error('Not found: %s', resultsRoot); end

    LABEL  = [ 25  65 125 165 285 325 365];
    ACTUAL = [365 325 285 165 125  65  25];

    % ── top-level README ─────────────────────────────────────────────────
    readmePath = fullfile(resultsRoot,'HEIGHT_LABELS_README.txt');
    L = strings(0,1);
    L(end+1) = "======================================================================";
    L(end+1) = " DROP-HEIGHT LABELS - READ BEFORE USING THIS DATA";
    L(end+1) = "======================================================================";
    L(end+1) = "";
    L(end+1) = "GB/shallow drop-height labels were recorded in REVERSE ORDER.";
    L(end+1) = "";
    L(end+1) = "Filenames, folder names, trialTags and meta.dropHeight_mm still carry";
    L(end+1) = "the ORIGINAL (incorrect) label. Nothing has been renamed, on purpose:";
    L(end+1) = "renaming would break trialTag references and, because the map is a";
    L(end+1) = "self-inverse swap, would leave no way to tell corrected from raw.";
    L(end+1) = "";
    L(end+1) = "CORRECTION MAP  (GB/shallow ONLY):";
    L(end+1) = "    labelled mm  ->  ACTUAL mm";
    for i = 1:numel(LABEL)
        L(end+1) = sprintf("    %8d     ->  %6d", LABEL(i), ACTUAL(i)); %#ok<AGROW>
    end
    L(end+1) = "";
    L(end+1) = "NOT affected: GB/full, CHIN/as_poured, CHIN/dense.";
    L(end+1) = "";
    L(end+1) = "EVIDENCE:";
    L(end+1) = "  With labelled heights, GB/shallow v0 / sqrt(2gh) reached 3.9-4.2.";
    L(end+1) = "  A ratio above 1 is impossible (faster than free fall). Measured v0";
    L(end+1) = "  was also nearly constant (~275, 265, 240 cm/s) across a 5x range of";
    L(end+1) = "  labelled heights. After correction all ratios fall to ~1.0-1.1.";
    L(end+1) = "";
    L(end+1) = "HOW TO USE THE DATA CORRECTLY:";
    L(end+1) = "  src/true_drop_height.m    single source of truth for the map";
    L(end+1) = "  src/load_meta_checked.m   use INSTEAD of load(tp,'meta'); it adds";
    L(end+1) = "                            meta.dropHeight_true_mm and warns once";
    L(end+1) = "";
    L(end+1) = "  Use meta.dropHeight_true_mm for ALL physics (v0, scaling, plots).";
    L(end+1) = "  Use meta.dropHeight_mm / trialTag ONLY as file identifiers.";
    L(end+1) = "";
    L(end+1) = "Depth (d_final_cm) is NOT affected by this correction.";
    L(end+1) = "The 4.00 cm model-overdepth and 2.50 cm shallow-bed rules are unchanged.";
    L(end+1) = "";
    L(end+1) = sprintf("written %s by scripts/write_height_labels.m", datestr(now));
    L(end+1) = "======================================================================";
    readmeTxt = strjoin(L, newline);

    % ── per-folder markers: any container folder holding shallow tracks ──
    D = dir(fullfile(resultsRoot,'**','*_shallow_tracks.mat'));
    folders = unique(string({D.folder}));
    folders = arrayfun(@(f) string(fileparts(f)), folders);   % tracks/ -> container/
    folders = cellstr(unique(folders(:)));   % cellstr: fopen needs char, not string

    markerTxt = strjoin([ ...
        "!! GB/shallow DROP-HEIGHT LABELS IN THIS FOLDER ARE REVERSED !!"; ""; ...
        "The height in the folder name / filenames / trialTags is the ORIGINAL"; ...
        "label and is NOT the physical drop height."; ""; ...
        "    labelled  25 ->  365 mm        labelled 285 ->  125 mm"; ...
        "    labelled  65 ->  325 mm        labelled 325 ->   65 mm"; ...
        "    labelled 125 ->  285 mm        labelled 365 ->   25 mm"; ...
        "    labelled 165 ->  165 mm  (unchanged)"; ""; ...
        "Use src/load_meta_checked.m, which supplies meta.dropHeight_true_mm."; ...
        "Full explanation: 03_RESULTS/HEIGHT_LABELS_README.txt"; ...
        "Depth values are unaffected."], newline);

    % ── act ──────────────────────────────────────────────────────────────
    if o.remove
        n = 0;
        if isfile(readmePath), delete(readmePath); n = n+1; end
        for k = 1:numel(folders)
            mp = fullfile(folders{k},'_HEIGHT_REVERSED.txt');
            if isfile(mp), delete(mp); n = n+1; end
        end
        fprintf('Removed %d marker file(s).\n', n);
        return;
    end

    fprintf('README : %s\n', readmePath);
    fprintf('markers: %d shallow container folder(s)\n', numel(folders));
    for k = 1:numel(folders)
        fprintf('   %s\n', fullfile(folders{k},'_HEIGHT_REVERSED.txt'));
    end
    if o.dryRun
        fprintf('\nDRY RUN - nothing written.\n');
        return;
    end

    fid = fopen(char(readmePath),'w');
    if fid < 0, error('Could not write %s', readmePath); end
    fprintf(fid,'%s\n',readmeTxt); fclose(fid);
    for k = 1:numel(folders)
        mp  = fullfile(folders{k},'_HEIGHT_REVERSED.txt');
        fid = fopen(mp,'w');
        if fid < 0, warning('Could not write %s', mp); continue; end
        fprintf(fid,'%s\n',markerTxt); fclose(fid);
    end
    fprintf('\nWrote README + %d marker file(s). No data touched.\n', numel(folders));
end
