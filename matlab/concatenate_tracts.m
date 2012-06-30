function concat_tracts = concatenate_tracts(tracts)

  num_tracts = size(tracts,1);
  
  concat_tracts = [];
  
  for tract_i = 1:num_tracts
    concat_tracts = [concat_tracts; tracts{tract_i,1}];
    concat_tracts = [concat_tracts; tracts{tract_i,2}];
    concat_tracts = [concat_tracts; tracts{tract_i,3}];    
  end

end