function trials = load_selected_trials()
% LOAD_SELECTED_TRIALS  Prompt user to select kinematics .mat files,
%   then optionally remove outlier trials before returning.

    trials = struct([]);
    i      = 0;

    fprintf('Select kinematics .mat files. Cancel dialog when done.\n');

    while true
        [file, path] = uigetfile('*_kinematics.mat', ...
                                  sprintf('Select trial %d (cancel when done)', i+1));
        if isequal(file, 0), break; end

        i = i + 1;
        S = load(fullfile(path, file));
        trials(i).trialInfo  = S.trialInfo;
        trials(i).kinematics = S.kinematics;
        trials(i).scalars    = S.scalars;
        trials(i).h_cm       = S.trialInfo.h_cm;
        trials(i).label      = sprintf('h = %.2f cm', S.trialInfo.h_cm);
        fprintf('Loaded: %s  (h=%.2f cm, v0=%.2f cm/s, d=%.3f cm)\n', ...
            file, S.trialInfo.h_cm, S.scalars.v0_cm_s, S.scalars.d_final_cm);
    end

    if isempty(trials), error('No files selected.'); end

    % Sort by h_cm
    [~, idx] = sort([trials.h_cm]);
    trials   = trials(idx);
    fprintf('Loaded %d trials across %d drop heights.\n', ...
        numel(trials), numel(unique([trials.h_cm])));

    % Optional removal 
    fprintf('\n── Loaded trials ────────────────────────────────────\n');
    for k = 1:numel(trials)
        fprintf('  [%2d]  h=%5.2f cm  |  T%02d  |  v0=%6.1f cm/s  |  d=%.3f cm\n', ...
            k, trials(k).h_cm, trials(k).trialInfo.trialNum, ...
            trials(k).scalars.v0_cm_s, trials(k).scalars.d_final_cm);
    end

    remove = input('\nEnter trial indices to remove (e.g. [3 7]), or [] to keep all: ');
    if ~isempty(remove)
        remove = unique(round(remove));
        remove = remove(remove >= 1 & remove <= numel(trials));
        trials(remove) = [];
        fprintf('Removed %d trial(s). %d remaining.\n', numel(remove), numel(trials));
    else
        fprintf('No trials removed.\n');
    end
end