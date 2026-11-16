function figs = paperFigures()
%PAPERFIGURES Names of the paper-result tests and their compute functions.
names = {'Fig1', 'Fig2', 'Fig3', 'Fig4_Fig5', 'Fig6a', 'Fig6b', 'Fig7', 'Fig8', 'Fig9_Fig10'};
figs = struct('name', names, 'compute', strcat('compute', names));
end
