function [tract_matrix, labels] = load_unzip_tracts(filename)
  
  [tracts, props, prop_keys, prop_values] = load_tracts(filename);
  [tract_matrix, labels] = unzip_tracts(tracts, prop_keys, prop_values);
  
end