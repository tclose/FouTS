function add_colour_key(bundle_indices, colours)

global colours_of_bundles

if ~exist('bundle_indices','var')
  bundle_indices=[0:size(colours_of_bundles,1)-1];
end

if isempty(bundle_indices)
  return
end

if exist('colours','var')
  
  if isempty(colours) || (max(bundle_indices) >= size(colours,1))
    set_bundle_colours(max(bundle_indices));  
    colours = colours_of_bundles;
  end
  
else  
  colours = colours_of_bundles;
end


unique_indices = sort(unique(bundle_indices));

h = figure();

set(h,'Units','normalized');
set(h, 'Position', [0.0 0.05 0.05 0.9]);
set(h, 'Name', 'Colour Key');

colours2 = reshape(colours, size(colours,1),1,3);

image(1, [0:1:(length(unique_indices)-1)], colours2((unique_indices+1),:,:));

set(gca, 'YTick', [0:1:(length(unique_indices)-1)]);
set(gca, 'YTickLabel', unique_indices);
set(gca, 'XTick', []);
set(get(gca,'YLabel'),'String','Bundle Index');