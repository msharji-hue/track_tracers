function sha = pipeline_commit(codeDir)
%PIPELINE_COMMIT  Short git commit of the pipeline, or '' if unavailable.
%
%   Frames are no longer kept, so a trial's detections are reproducible only
%   from raw video + code + parameters. This records which code.
%
%   Cheap by construction: the answer is cached in a persistent variable, so
%   the shell is invoked at most once per MATLAB session however many trials a
%   batch processes. Any failure -- no git on PATH, not a repository, a
%   detached or shallow checkout -- returns '' rather than raising, because
%   provenance is a nice-to-have and must never cost a trial.
%
%   A '-dirty' suffix means the working tree had uncommitted changes when the
%   trial was processed, so the commit alone does not identify the code that
%   ran. Treat those results as unversioned.

    persistent cached
    if ~isempty(cached), sha = cached; return; end

    sha = '';
    try
        if nargin < 1 || isempty(codeDir)
            codeDir = fileparts(fileparts(mfilename('fullpath')));
        end
        cmd = sprintf('git -C "%s" rev-parse --short HEAD', codeDir);
        [st, out] = system(cmd);
        if st == 0
            sha = strtrim(out);
            % Uncommitted changes mean the sha does not pin what ran.
            [st2, out2] = system(sprintf('git -C "%s" status --porcelain', codeDir));
            if st2 == 0 && ~isempty(strtrim(out2))
                sha = [sha '-dirty'];
            end
        end
    catch
        sha = '';
    end
    cached = sha;
end
