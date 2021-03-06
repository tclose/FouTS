function colour_indices = add_colour_key(bundle_indices, colours, colour_indices)

global colours_of_bundles

if ~exist('bundle_indices', 'var')
    bundle_indices = 0:1:(length(colours_of_bundles)-1);
end

if ~exist('colours', 'var')
    colours = colours_of_bundles;
end

if ~exist('colour_indices', 'var')
    [~, colour_indices] = set_bundle_colours(bundle_indices, colours);
end
    
unique_indices = sort(unique(bundle_indices));

h = figure();

set(h,'Units','normalized');
set(h, 'Position', [0.0 0.05 0.05 0.9]);
set(h, 'Name', 'Colour Key');

colours2 = reshape(colours, size(colours,1),1,3);

image(1, 0:1:(length(unique_indices)-1), colours2((colour_indices(unique_indices + 1)),:,:));

set(gca, 'YTick', [0:1:(length(unique_indices)-1)]);
set(gca, 'YTickLabel', unique_indices);
set(gca, 'XTick', []);
set(get(gca,'YLabel'),'String','Bundle Index');