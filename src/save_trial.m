function save_trial(trialInfo, kinematics, scalars, figs, figNames, outRoot)
% SAVE_TRIAL  Print summary, save kinematics .mat/.csv, and export figures.
%
%   Inputs:
%       trialInfo  - struct: material, batchName, heightLabel, trialNum, h_cm
%       kinematics - struct: full kinematics + geometry fields
%       scalars    - struct: v0_cm_s, d_final_cm, t_stop_s, h_cm, H_cm
%       figs       - array of figure handles to save (e.g. [1 2 3 4])
%       figNames   - cell array of base name strings, one per figure handle
%                    e.g. {'kinematics','tracking_qa','phase_space','impact_qa'}
%       outRoot    - results root folder path

    % ── Summary ───────────────────────────────────────────────────────────
    fprintf('\n── Trial summary ────────────────────────────────\n');
    fprintf('  %s  T%02d\n', trialInfo.heightLabel, trialInfo.trialNum);
    fprintf('  v0       = %.2f cm/s\n',  scalars.v0_cm_s);
    fprintf('  d_final  = %.3f cm\n',    scalars.d_final_cm);
    fprintf('  H        = %.3f cm\n',    scalars.H_cm);
    fprintf('  t_stop   = %.4f s\n',     scalars.t_stop_s);
    fprintf('─────────────────────────────────────────────────\n\n');

    % ── Kinematics .mat + scalars .csv ────────────────────────────────────
    save_kinematics(trialInfo, kinematics, scalars, outRoot);

    % ── Figures ───────────────────────────────────────────────────────────
    figDir = fullfile(outRoot, 'figures', trialInfo.heightLabel);
    if ~exist(figDir, 'dir'), mkdir(figDir); end

    for k = 1:numel(figs)
        fh      = figure(figs(k));
        figBase = sprintf('%s_T%02d_%s', ...
            trialInfo.heightLabel, trialInfo.trialNum, figNames{k});
        exportgraphics(fh, fullfile(figDir, [figBase '.png']), 'Resolution', 300);
        savefig(fh,         fullfile(figDir, [figBase '.fig']));
        fprintf('Saved figure: %s\n', figBase);
    end
end
