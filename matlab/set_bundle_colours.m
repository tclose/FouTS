function set_bundle_colours(max_bundle_index)

  global colours_of_bundles;

  num_colours = max_bundle_index+1;
  
  if size(colours_of_bundles,1) < num_colours
    
		colours_of_bundles = rand([num_colours, 3]);
% 		add_colour_key(colours_of_bundles, 0:1:num_colours-1);
    
  end

end