function tracts = parse_tracts(elems)
  
  tract_divs = split_at_file_seperator(elems, [-inf,-inf,-inf]);
  
  if isempty(tract_divs)
    tract_divs = split_at_file_seperator(elems, [-inf, -inf, inf]); % Legacy issue with incorrect seperators.
  end
  
  tracts = cell(0);
  
  for tract_i = 1:size(tract_divs,1)
    
    tract_row= split_at_file_seperator(tract_divs{tract_i}, [nan,nan,nan]);
    
    if size(tract_row,1) ~= 3
      error([num2str(size(tract_row,1)) ' strands were loaded from ' num2str(tract_i) ' instead of 3']);
    end
    
    for ax_i = 1:3
      tracts{tract_i,ax_i} = tract_row{ax_i}; 
    end
    
  end
  
end

