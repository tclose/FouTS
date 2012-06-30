function [tcks, props, prop_keys, prop_values] = load_tcks(filename)

  [elems, props] = read_fibres(filename);

  tcks = parse_strands(elems);
  
  if exist([filename 'x']) == 2
    [prop_keys, prop_values] = read_element_properties([filename 'x']);
  else
    prop_keys = cell(0);
    prop_values = cell(0);
  end
  
  
  if (size(prop_values,1) ~= size(tcks,1))  && size(prop_values,1) ~= 0
    error (['Number of loaded property values (' num2str(size(prop_values,1)) ') does not match number of loaded tracks (' num2str(size(tcks,1)) ').']);
  end
    
  
end