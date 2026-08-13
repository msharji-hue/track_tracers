function inspect_impact(root, varargin)
%INSPECT_IMPACT  Plot z(t) around the assigned impact frame.
%
%   Diagnostic only: reads *_kin.mat, writes nothing, changes nothing.
%
%   The question this answers: at the assigned impact index, is there a
%   pre-impact free-fall segment in the record, or does the trace already
%   start mid-event? If impact sits at frame ~3 of the usable record, the
%   depth d = z - z(impact) is truncated and the trial cannot yield a valid
%   penetration depth.
%
%   What to look for in each panel:
%     - GREY  region left of the dashed line = frames before assigned impact.
%       Wide grey band with a straight, steepening fall  -> impact captured.
%       No grey band, trace begins already curving       -> impact missed.
%     - The red dotted line is free fall from the true drop height. The
%       pre-impact segment should run parallel to it.
%
%   USAGE
%       root = 'D:\ME_GRANULAB\JerboaImpact';
%       inspect_impact(root);                       % default: the suspect set
%       inspect_impact(root, 'Tags', ["165mm_T01_dense","365mm_T04_dense"]);
%
%   OPTIONS
%       'Tags'     explicit trial tags to plot (overrides the default set)
%       'PerCase'  how many trials per default case (default 2)

opt.Tags    = strings(0);
opt.PerCase = 2;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

G = 980;

% ---------------------------------------------------------------- locate
F = dir(fullfile(root, '03_RESULTS', '**', '*_kin.mat'));
F = F(~[F.isdir]);
if isempty(F)
    error('inspect_impact:noFiles', 'No *_kin.mat under %s', ...
          fullfile(root,'03_RESULTS'));
end
allTags = string(erase({F.name}, '_kin.mat'))';

if isempty(opt.Tags)
    % two suspect CHIN/dense cases + one GB/full control at the same height
    want = {"165mm_", "_dense",  'CHIN/dense  h = 16.5 cm  (suspect)'; ...
            "365mm_", "_dense",  'CHIN/dense  h = 36.5 cm  (suspect)'; ...
            "365mm_", "_full",   'GB/full     h = 36.5 cm  (control)'};
    tags = strings(0); labs = strings(0);
    for k = 1:size(want,1)
        hit = allTags(startsWith(allTags, want{k,1}) & endsWith(allTags, want{k,2}));
        hit = hit(1:min(opt.PerCase, numel(hit)));
        tags = [tags; hit];                       %#ok<AGROW>
        labs = [labs; repmat(string(want{k,3}), numel(hit), 1)]; %#ok<AGROW>
    end
else
    tags = string(opt.Tags(:));
    labs = strings(numel(tags),1);
end
if isempty(tags)
    error('inspect_impact:noMatch', 'No matching trials found.');
end

% ---------------------------------------------------------------- plot
n  = numel(tags);
nc = min(3, n);
nr = ceil(n/nc);
figure('Color','w','Position',[60 60 420*nc 340*nr]);
tl = tiledlayout(nr, nc, 'Padding','compact', 'TileSpacing','compact');

for i = 1:n
    k  = find(allTags == tags(i), 1);
    S  = load(fullfile(F(k).folder, F(k).name));
    kin = S.kin;  meta = S.meta;

    if i == 1
        fprintf('kin fields : %s\n',  strjoin(fieldnames(kin)',  ', '));
        fprintf('meta fields: %s\n\n', strjoin(fieldnames(meta)', ', '));
    end

    z   = local_pick(kin, {'z_cm','z','depth_cm','zRod_cm','pos_cm'});
    fps = local_pick(meta, {'fps_true','fps'});
    if isempty(fps), fps = local_pick(kin, {'fps'}); end
    iImp = local_pick(kin, {'impact_frame','impactIdx','iImpact','impact_index'});
    iStp = local_pick(kin, {'stop_frame','stopIdx','iStop','stop_index'});
    if isempty(z) || isempty(iImp)
        error('inspect_impact:fields', ...
             ['Could not find z or impact index in %s.\n' ...
              'kin fields are: %s'], tags(i), strjoin(fieldnames(kin)', ', '));
    end
    if isempty(fps), fps = 2778; end

    z = z(:);  t = (0:numel(z)-1)'/fps * 1000;    % ms
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    % shade the pre-impact record
    if iImp > 1
        patch(ax, [t(1) t(iImp) t(iImp) t(1)], ...
              [min(z) min(z) max(z) max(z)], [0.90 0.90 0.90], ...
              'EdgeColor','none', 'HandleVisibility','off');
    end
    plot(ax, t, z, '-', 'Color',[0.1 0.3 0.7], 'LineWidth',1.4);
    xline(ax, t(iImp), 'k--', sprintf('impact (frame %d)', iImp), ...
          'LabelOrientation','horizontal', 'FontSize',8);
    if ~isempty(iStp) && iStp <= numel(t)
        xline(ax, t(iStp), 'r:', sprintf('stop (frame %d)', iStp), ...
              'LabelOrientation','horizontal', 'FontSize',8);
    end

    % free-fall reference through the impact point
    hmm = local_pick(meta, {'dropHeight_true_mm','dropHeight_mm'});
    if ~isempty(hmm)
        vff = sqrt(2*G*hmm/10);
        tt  = t(max(1,iImp-25):iImp)/1000;
        plot(ax, tt*1000, z(iImp) - vff*(tt - t(iImp)/1000), 'r:', 'LineWidth',1.2);
    end

    ttl = tags(i); if labs(i) ~= "", ttl = labs(i) + newline + tags(i); end
    title(ax, ttl, 'FontSize',9, 'Interpreter','none');
    xlabel(ax,'t (ms)'); ylabel(ax,'z (cm)');

    fprintf('%-26s impact frame %4d of %5d   pre-impact frames: %4d %s\n', ...
        tags(i), iImp, numel(z), iImp-1, ...
        local_tern(iImp-1 < 10, '  <-- TOO FEW, trace starts mid-event', ''));
end

title(tl, ['z(t) around the assigned impact.  Grey = record before impact.  ' ...
           'Red dotted = free fall from true height.'], 'FontWeight','bold','FontSize',9);
end

% ------------------------------------------------------------------ helpers
function v = local_pick(S, names)
v = [];
for i = 1:numel(names)
    if isfield(S, names{i}), v = S.(names{i}); return; end
end
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
