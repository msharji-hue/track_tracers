function save_all_figures(figs, names, outDir)
% SAVE_ALL_FIGURES  Export figure handles to png, pdf, and fig subfolders.
%
%   Inputs:
%       figs   - cell array of figure handles
%       names  - cell array of base filenames (no extension)
%       outDir - root output directory

    pngDir = fullfile(outDir, 'png');
    pdfDir = fullfile(outDir, 'pdf');
    figDir = fullfile(outDir, 'fig');
    mkdir(pngDir); mkdir(pdfDir); mkdir(figDir);

    for f = 1:numel(figs)
        base = names{f};
        hFig = figs{f};

        if ~isgraphics(hFig, 'figure')
            warning('save_all_figures: invalid handle for "%s" — skipping.', base);
            continue;
        end

        figure(hFig); drawnow;

        try
            exportgraphics(hFig, fullfile(pngDir, [base '.png']), ...
                'Resolution', 300, 'BackgroundColor', 'white');
            exportgraphics(hFig, fullfile(pdfDir, [base '.pdf']), ...
                'ContentType', 'vector', 'BackgroundColor', 'white');
            savefig(hFig, fullfile(figDir, [base '.fig']));
            fprintf('Saved: %s\n', base);
        catch ME
            warning('save_all_figures: failed to save "%s": %s', base, ME.message);
        end
    end

    fprintf('\nAll figures saved to:\n  %s\n', outDir);
end
