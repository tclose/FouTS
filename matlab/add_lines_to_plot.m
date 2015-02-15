function add_lines_to_plot(tcks, colours, bundle_indices, line_style)
% Adds lines to voxel plot.

 if ~exist('line_style', 'var')
     line_style = '-';
 end

  num_tcks = length(tcks);

  hold on;
  
  for tck_i = 1:num_tcks

    if ~isempty(tcks{tck_i})
        h = plot3(tcks{tck_i}(:,1), tcks{tck_i}(:,2), tcks{tck_i}(:,3), line_style);

        set(h, 'Color', colours(bundle_indices(tck_i)+1,:));
    end
  end
  
  hold off;

end
