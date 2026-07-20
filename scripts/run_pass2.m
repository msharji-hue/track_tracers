function run_pass2(retryListFile, varargin)
% RUN_PASS2  Re-run ONLY the no-8-marker failures with a fixed detection
% setting, without disturbing the 341 already-tracked trials.
%
%   Forces sensitivity and edgeThresh to the Pass 2 values, then calls the
%   normal batch on the supplied failure list with policy='overwrite' (so the
%   previously-failed trials are re-attempted; already-successful trials are
%   NOT in the list and are never touched). The 8-marker rule is unchanged:
%   a trial is recovered only if a frame with EXACTLY 8 markers is found;
%   otherwise it simply fails again and is left as-is.
%
%   Frames go to JERBOA_FRAMES_ROOT scratch (set in the shell before calling),
%   so no frame PNGs land in Dropbox.
%
%   run_pass2(retryListFile)
%   run_pass2(retryListFile, 'sensitivity',0.87, 'edgeThresh',0.10, ...
%                            'outputRoot',OUT, 'batchLabel','Batch 5')
%
%   retryListFile : text file with one absolute .avi path per line (the
%                   no-8-marker failures). See the awk command in the notes.

    p = inputParser;
    addParameter(p,'sensitivity',0.87);
    addParameter(p,'edgeThresh', 0.10);
    addParameter(p,'outputRoot', '/Users/muhannadalsharji/Library/CloudStorage/Dropbox-UniversityofMichigan/Muhannad Al Sharji/ME_GRANULAB/JerboaImpact');
    addParameter(p,'batchLabel', 'Batch 5');
    parse(p,varargin{:});
    o = p.Results;

    if ~isfile(retryListFile)
        error('run_pass2: list file not found: %s', retryListFile);
    end

    fprintf('=== PASS 2 (fixed detection settings) ===\n');
    fprintf('  list        : %s\n', retryListFile);
    fprintf('  sensitivity : %.2f  (default is 0.85)\n', o.sensitivity);
    fprintf('  edgeThresh  : %.2f\n', o.edgeThresh);
    fprintf('  policy      : overwrite (only the listed failures)\n');
    fprintf('  8-marker rule: UNCHANGED (exactly 8 required)\n\n');

    % Build options for the standard batch, injecting the params override.
    opts = struct( ...
        'videoListFile', retryListFile, ...
        'policy',        'overwrite', ...
        'batchLabel',    o.batchLabel, ...
        'detectParams',  struct('sensitivity',o.sensitivity,'edgeThresh',o.edgeThresh));

    process_trial('batch', '', o.outputRoot, opts);
end
