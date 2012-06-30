function [strand_matrix, labels] = load_unzip_strands(filename)
  
  [strands, props, prop_keys, prop_values] = load_strands(filename);
  [strand_matrix, labels] = unzip_strands(strands, prop_keys, prop_values);
  
end