function save_analysis(heights, hVals, cmap, outDir)
% SAVE_ANALYSIS  Save grouped trial data for downstream analysis.
%   Saves analysis_data.mat to outDir.
%
%   Inputs:
%       heights - struct array from group_trials_by_height
%       hVals   - array of unique h_cm values
%       cmap    - [nH x 3] colormap used in overlay_plots
%       outDir  - folder to save analysis_data.mat

    v0_all    = [];
    d_all     = [];
    tstop_all = [];
    h_all     = [];
    agmax_all = [];

    for j = 1:numel(heights)
        hg = heights(j);
        v0_all    = [v0_all,    arrayfun(@(t) t.scalars.v0_cm_s,    hg.trials)];
        d_all     = [d_all,     arrayfun(@(t) t.scalars.d_final_cm, hg.trials)];
        tstop_all = [tstop_all, arrayfun(@(t) t.scalars.t_stop_s,   hg.trials)];
        h_all     = [h_all,     repmat(hg.h_cm, 1, hg.nTrials)];
        agmax_all = [agmax_all, arrayfun(@(t) max( ...
            t.kinematics.a_plus_g(t.kinematics.impact_index:t.kinematics.stopFrame), ...
            [], 'omitnan'), hg.trials)];
    end

    analysisData = struct( ...
        'heights',   heights,   ...
        'hVals',     hVals,     ...
        'cmap',      cmap,      ...
        'v0_all',    v0_all,    ...
        'd_all',     d_all,     ...
        'tstop_all', tstop_all, ...
        'h_all',     h_all,     ...
        'agmax_all', agmax_all);

    outMat = fullfile(outDir, 'analysis_data.mat');
    save(outMat, 'analysisData', '-v7.3');
    fprintf('Saved analysis_data.mat: %d trials across %d heights\n', ...
        numel(v0_all), numel(heights));
end
