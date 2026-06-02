location = ['.\Figures'];

name = "Fig6a";

locname = fullfile(location, name);

savefig(locname);
exportgraphics(gcf, strcat(locname, ".png"));
exportgraphics(gcf, strcat(locname, ".eps"));
exportgraphics(gcf, strcat(locname, ".emf"));
exportgraphics(gcf, strcat(locname, ".pdf"));

