function [tracts, props, prop_keys, prop_values] = load_tracts(filename)

  [elems, props] = read_fibres(filename);

  tracts = parse_tracts(elems);

  if exist([filename 'x']) == 2
    [prop_keys, prop_values] = read_element_properties([filename 'x']);

    if size(prop_values,1) ~= size(tracts,1)
      error (['Number of loaded property values (' num2str(size(prop_values,1)) ') does not match number of loaded tracts (' num2str(size(tracts,1)) ').']);
    end  

  else
    prop_keys = cell(0);
    prop_values = cell(0);
  end
  
end

