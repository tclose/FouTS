function [section_matrix, labels] = load_unzip_tract_sections(filename)

  [sections, prop_keys, prop_values] = load_tract_sections(filename);
  [section_matrix, labels]           = unzip_tract_sections(sections, prop_keys, prop_values);
  
end