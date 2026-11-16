function labels = concentrationLabels(species, c)
%CONCENTRATIONLABELS Legend entries '$c_\mathrm{<species>} = <c>$ [wppm]'.
labels = arrayfun(@(ci) sprintf('$c_\\mathrm{%s} = %d$ [wppm]', species, ci), ...
    c, 'UniformOutput', false);
end
