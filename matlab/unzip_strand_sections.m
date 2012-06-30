function [section_matrix, labels] = unzip_strand_sections(sections, prop_keys, prop_values)
  
  num_sections = size(sections,1);

  if num_sections == 0
    error ('No sections loaded.');
  end

  num_props = length(prop_keys);    
  num_parameters = 6 + num_props;

  section_matrix = zeros(num_parameters, num_sections);
  labels = cell(num_parameters,1);

  for prop_i = 1:num_props
    section_matrix(prop_i,:) = get_properties(prop_keys, prop_values, prop_keys{prop_i});
    labels{prop_i} = prop_keys{prop_i};
  end

  for section_i = 1:num_sections
    section_matrix((num_props+1):end, section_i) = sections{section_i}(:);
  end
  
  labels{num_props+1} = 'Position - Dim 0';
  labels{num_props+2} = 'Position - Dim 1';
  labels{num_props+3} = 'Position - Dim 2';
  labels{num_props+4} = 'Tangent - Dim 0';
  labels{num_props+5} = 'Tangent - Dim 1';
  labels{num_props+6} = 'Tangent - Dim 2';
  
end