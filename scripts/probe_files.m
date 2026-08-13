function probe_files(root)
%PROBE_FILES  Report what is actually inside the saved files. READ-ONLY.
%
%   Diagnostic for scan_impact_dist failing to match. Prints structure only;
%   loads at most two files and writes nothing.
%
%   USAGE
%       probe_files('D:\ME_GRANULAB\JerboaImpact')
%
%   Paste the whole output back and the scanner can be corrected exactly.

fprintf('\n================ probe_files ================\n');

% ---------------------------------------------------------------- 1. counts
TK = dir(fullfile(root,'03_RESULTS','**','*_tracks.mat')); TK = TK(~[TK.isdir]);
KN = dir(fullfile(root,'03_RESULTS','**','*_kin.mat'));    KN = KN(~[KN.isdir]);
CS = dir(fullfile(root,'03_RESULTS','**','*_kin_scalars.csv')); CS = CS(~[CS.isdir]);
fprintf('\n1. FILE COUNTS\n');
fprintf('   *_tracks.mat       : %d\n', numel(TK));
fprintf('   *_kin.mat          : %d\n', numel(KN));
fprintf('   *_kin_scalars.csv  : %d\n', numel(CS));
if isempty(TK)
    fprintf('\n   No tracks found under %s\n', fullfile(root,'03_RESULTS'));
    fprintf('   Listing two levels down so we can see the real layout:\n');
    d = dir(fullfile(root,'03_RESULTS','*')); d = d([d.isdir] & ~startsWith({d.name},'.'));
    for i = 1:min(6,numel(d))
        fprintf('     %s\n', d(i).name);
        e = dir(fullfile(d(i).folder,d(i).name,'*'));
        e = e([e.isdir] & ~startsWith({e.name},'.'));
        for j = 1:min(4,numel(e)), fprintf('        %s\n', e(j).name); end
    end
    return
end
fprintf('   example tracks path: %s\n', fullfile(TK(1).folder, TK(1).name));

% ------------------------------------------------- 2. inside a _tracks.mat
fprintf('\n2. INSIDE %s\n', TK(1).name);
S = load(fullfile(TK(1).folder, TK(1).name));
top = fieldnames(S);
fprintf('   top-level variables: %s\n', strjoin(top', ', '));
for i = 1:numel(top)
    v = S.(top{i});
    if isstruct(v)
        fprintf('   %s is a struct with fields:\n', top{i});
        f = fieldnames(v);
        for j = 1:numel(f)
            x = v.(f{j});
            fprintf('      %-22s %-10s %s\n', f{j}, class(x), mat2str(size(x)));
        end
    else
        fprintf('   %-22s %-10s %s\n', top{i}, class(v), mat2str(size(v)));
    end
end

% ---------------------------------------------------- 3. inside a _kin.mat
if ~isempty(KN)
    fprintf('\n3. INSIDE %s\n', KN(1).name);
    K = load(fullfile(KN(1).folder, KN(1).name));
    ktop = fieldnames(K);
    fprintf('   top-level variables: %s\n', strjoin(ktop', ', '));
    for i = 1:numel(ktop)
        v = K.(ktop{i});
        if ~isstruct(v), continue; end
        fprintf('   %s fields:\n', ktop{i});
        f = fieldnames(v);
        for j = 1:numel(f)
            x = v.(f{j});
            extra = '';
            if isnumeric(x) && isscalar(x), extra = sprintf('= %g', x); end
            if (ischar(x)||isstring(x)) && numel(string(x))==1, extra = sprintf('= %s', string(x)); end
            fprintf('      %-22s %-10s %-12s %s\n', f{j}, class(x), mat2str(size(x)), extra);
        end
    end
end

% ------------------------------------------------- 4. condition strings
fprintf('\n4. CONDITION VALUES IN _kin_scalars.csv\n');
cond = strings(numel(CS),1); tags = strings(numel(CS),1);
for i = 1:numel(CS)
    T = readtable(fullfile(CS(i).folder, CS(i).name));
    if height(T)<1, continue; end
    tags(i) = string(erase(CS(i).name,'_kin_scalars.csv'));
    if ismember('condition', T.Properties.VariableNames)
        cond(i) = strtrim(string(T.condition(1)));
    end
end
u = unique(cond(cond~="" ));
for i = 1:numel(u)
    fprintf('   "%s"   (n = %d)\n', u(i), sum(cond==u(i)));
end
if isempty(u), fprintf('   no condition column found\n'); end

% ------------------------------------------------- 5. tag matching
fprintf('\n5. TAG MATCHING  tracks <-> scalars\n');
tkTags = string(erase({TK.name},'_tracks.mat'))';
fprintf('   example tracks tag : "%s"\n', tkTags(1));
fprintf('   example scalars tag: "%s"\n', tags(find(tags~="",1)));
fprintf('   tracks tags that match a scalars tag: %d of %d\n', ...
    sum(ismember(tkTags, tags)), numel(tkTags));

fprintf('\n============ end probe_files ============\n\n');
end
