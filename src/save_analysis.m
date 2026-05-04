function save_analysis(heights, hVals, cmap, peakZ, peakV, peakAG, outDir)
% SAVE_ANALYSIS  Save all data needed for downstream normalization and fitting.
%   Saves analysis_data.mat to outDir.
%
%   Load in another script with:
%       [file, path] = uigetfile('analysis_data.mat', 'Select analysis data');
%       load(fullfile(path, file));
%       % then use analysisData.v0_all, analysisData.heights, etc.
%
%   Inputs:
%       heights  - struct array from group_trials_by_height
%       hVals    - array of unique h_cm values
%       cmap     - [nH x 3] colormap used in overlay_plots
%       peakZ    - d_final inset data from Fig 1
%       peakV    - v0 inset data from Fig 2
%       peakAG   - peak a+g inset data from Fig 3
%       outDir   - folder to save analysis_data.mat

    % ── Flat arrays across all trials ─────────────────────────────────────
    v0_all    = [];
    d_all     = [];
    tstop_all = [];
    h_all     = [];
    agmax_all = [];

    for j = 1:numel(heights)
        hg = heights(j);

        v0_all    = [v0_all,    arrayfun(@(t) t.scalars.v0_cm_s,    hg.trials)];
        d_all     = [d_all,     arrayfun(@(t) t.scalars.d_final_cm,  hg.trials)];
        tstop_all = [tstop_all, arrayfun(@(t) t.scalars.t_stop_s,    hg.trials)];
        h_all     = [h_all,     repmat(hg.h_cm, 1, hg.nTrials)];
        agmax_all = [agmax_all, arrayfun(@(t) max( ...
                        t.kinematics.a_plus_g( ...
                            t.kinematics.impact_index:t.kinematics.stopFrame), ...
                        [], 'omitnan'), hg.trials)];
    end

    % ── Pack into single struct ───────────────────────────────────────────
    analysisData = struct( ...
        'heights',   heights,   ...   % full group struct — all kinematics
        'hVals',     hVals,     ...   % unique drop heights
        'cmap',      cmap,      ...   % colormap for consistent plotting
        'peakZ',     peakZ,     ...   % d_final vs v0 (Fig 1 inset)
        'peakV',     peakV,     ...   % v0 vs h     (Fig 2 inset)
        'peakAG',    peakAG,    ...   % peak a+g vs v0 (Fig 3 inset)
        'v0_all',    v0_all,    ...   % impact velocity per trial [1 x N]
        'd_all',     d_all,     ...   % final depth per trial     [1 x N]
        'tstop_all', tstop_all, ...   % stop time per trial       [1 x N]
        'h_all',     h_all,     ...   % drop height per trial     [1 x N]
        'agmax_all', agmax_all);      % peak a+g per trial        [1 x N]

    % ── Save ──────────────────────────────────────────────────────────────
    outMat = fullfile(outDir, 'analysis_data.mat');
    save(outMat, 'analysisData', '-v7.3');
    fprintf('Saved: %s\n  %d trials across %d heights\n', ...
        outMat, numel(v0_all), numel(heights));
end
