function tcks = strands2tcks(strands, num_length_sections)

  if ~exist('num_length_sections','var')
    num_length_sections = 100;
  end

  tcks = cell(size(strands));

  for strand_i = 1:length(strands) 
    tcks{strand_i} = fourier2tck(strands{strand_i},num_length_sections);
  end
  
end