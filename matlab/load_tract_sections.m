function [sections, prop_keys, prop_values] = load_tract_sections(filename)

  sections = parse_tract_sections(read_fibres(filename));
  
  
  [prop_keys, prop_values] = read_element_properties([filename 'x']);
  
  if size(prop_values,1) ~= size(sections,1)
    error (['Number of loaded property values (' num2str(size(prop_values,1)) ') does not match number of loaded segments (' num2str(size(sections,1)) ').']);
  end
    
  
end