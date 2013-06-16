function add_lines_to_plot(tcks, colours, bundle_indices)
% Adds lines to voxel plot.

  num_tcks = length(tcks);

  hold on;
  
  for tck_i = 1:num_tcks

    if ~isempty(tcks{tck_i})
        h = plot3(tcks{tck_i}(:,1), tcks{tck_i}(:,2), tcks{tck_i}(:,3));

        set(h, 'Color', colours(bundle_indices(tck_i)+1,:));
    end
  end
  
  hold off;

end
