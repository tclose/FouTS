function [num_length_sections, num_width_sections] = get_num_section_properties(properties, style)

  if isfield(properties, 'exp_num_segments')
    num_length_sections = str2double(properties.exp_num_segments);
  elseif isfield(properties, 'exp_num_length_sections')
    num_length_sections = str2double(properties.exp_num_length_sections);
  else
    disp('Warning!! Did not find number of length samples, using defaults.');
    num_length_sections = 50;
  end

  need_width_sections = false;

  if (strfind('lines', style) == 1)
    need_width_sections = true;
  elseif (strfind('tubes', style) == 1)
    need_width_sections = true;
  end

  if need_width_sections
    if isfield(properties, 'exp_num_strands')
      num_width_sections = str2double(properties.exp_num_strands);
    elseif isfield(properties, 'exp_num_width_sections')
      num_width_sections = str2double(properties.exp_num_width_sections);
    else
      disp('Warning!! Did not find number of width samples, using defaults.');
      num_width_sections = 2;
    end
  else
    num_width_sections = [];
  end
    
end
