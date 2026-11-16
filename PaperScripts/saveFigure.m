function saveFigure(name, fig)
%SAVEFIGURE Save a figure to PaperScripts/Figures as .fig, .png, .eps, .emf and .pdf.
%
%   saveFigure('Fig6a')       saves the current figure
%   saveFigure('Fig6a', fig)  saves figure fig

if nargin < 2
    fig = gcf;
end
location = fullfile(fileparts(mfilename('fullpath')), 'Figures', name);
savefig(fig, location);
for ext = {'.png', '.eps', '.emf', '.pdf'}
    exportgraphics(fig, [location ext{1}]);
end
end
