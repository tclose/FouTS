function [tract_sets_matrix, labels] = load_unzip_tract_sets(filename)
  
  [tract_sets, ~, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values] = load_tract_sets(filename);
  [tract_sets_matrix, labels] = unzip_tract_sets(tract_sets, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values);
  
end