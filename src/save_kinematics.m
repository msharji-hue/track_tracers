function save_kinematics(trialInfo, kinematics, scalars, outRoot)
% SAVE_KINEMATICS  Save all trial kinematics and publishable quantities.
%
%   Inputs:
%       trialInfo  - struct with material, batchName, heightLabel, trialNum, h_cm
%       kinematics - struct with full time series (t_s, depthRod_cm, v_smooth, etc.)
%       scalars    - struct with v0_cm_s, d_final_cm, t_stop_s
%       outRoot    - path to results root folder

    heightLabel = trialInfo.heightLabel;
    trialID     = sprintf('T%02d', trialInfo.trialNum);
    baseName    = sprintf('%s_%s', heightLabel, trialID);
    kinDir      = fullfile(outRoot, 'kinematics', heightLabel);
    if ~exist(kinDir, 'dir'), mkdir(kinDir); end

    % Save full kinematics .mat
    outMat = fullfile(kinDir, [baseName '_kinematics.mat']);
    save(outMat, 'trialInfo', 'kinematics', 'scalars', '-v7.3');
    fprintf('Saved: %s\n', outMat);

    % Save scalars to CSV for easy overlay/analysis
    outCsv = fullfile(kinDir, [baseName '_scalars.csv']);
    fid    = fopen(outCsv, 'w');
    fprintf(fid, 'material,batch,heightLabel,trialNum,h_cm,H_cm,v0_cm_s,d_final_cm,t_stop_s\n');
    fprintf(fid, '%s,%s,%s,%d,%.2f,%.2f,%.4f,%.4f,%.4f\n', ...
    trialInfo.material, trialInfo.batchName, trialInfo.heightLabel, ...
    trialInfo.trialNum, trialInfo.h_cm, scalars.H_cm, ...
    scalars.v0_cm_s, scalars.d_final_cm, scalars.t_stop_s);
    fclose(fid);
    fprintf('Saved: %s\n', outCsv);
end
