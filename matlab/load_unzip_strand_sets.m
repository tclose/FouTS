function [strand_sets_matrix, labels] = load_unzip_strand_sets(filename)
  
  [strand_sets, ~, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values] = load_strand_sets(filename);
  [strand_sets_matrix, labels] = unzip_strand_sets(strand_sets, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values);
  
end