function [tract_sets, properties, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values]  = load_tract_sets(filename, set_indices)

  [elems, properties] = read_fibres(filename);
  
  sets = split_at_file_seperator(elems, [-inf,nan,inf]);

  num_sets = length(sets);
  
  tract_sets = cell(num_sets,1);
  
  for set_i =1:num_sets
    tract_sets{set_i} = parse_tracts(sets{set_i});
  end
  
  if exist([filename 'x'],'file')
    [set_prop_keys, set_prop_values] = read_element_properties([filename 'x']);
  else
    set_prop_keys = [];
    set_prop_values = [];
  end
  
  if exist([filename 'xx'], 'file')
    [elem_prop_keys, elem_prop_values] = read_set_element_properties([filename 'xx']);  
  else
    elem_prop_keys = [];
    elem_prop_values = cell(num_sets,1);
  end
  
end