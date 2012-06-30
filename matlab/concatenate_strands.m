function concat_strands = concatenate_strands(strands)

  num_strands = size(strands,1);
  
  concat_strands = [];
  
  for strand_i = 1:num_strands
    
    concat_strands = [concat_strands; strands{strand_i}];
    
  end

end