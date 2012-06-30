function add_strands_to_plot(strands, radii, colours, bundle_indices, varargin)

  tcks = strands2tcks(strands);
  
  add_tcks_to_plot(tcks, radii, colours, bundle_indices, varargin{:}); 

end