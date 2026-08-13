function outRoot = resolve_output_root(preset)
% RESOLVE_OUTPUT_ROOT  Shared output-root resolver for the batch stages.
%   Returns preset if given; otherwise prompts once. '' on cancel.
%   (Extracted verbatim from process_trial.m so process_trial and
%   track_tracers_2 resolve the root identically.)
    if ~isempty(preset), outRoot = preset; return; end
    outRoot = uigetdir(pwd, ...
        'Select OUTPUT root (parent of 01_FRAMES / 02_SAVED_DETECTIONS / 03_RESULTS)');
    if isequal(outRoot, 0), outRoot = ''; end
end
