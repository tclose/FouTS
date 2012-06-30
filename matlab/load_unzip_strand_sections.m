function [section_matrix, labels] = load_unzip_strand_sections(filename)

  [sections, prop_keys, prop_values] = load_strand_sections(filename);
  [section_matrix, labels]           = unzip_strand_sections(sections, prop_keys, prop_values);
  
end