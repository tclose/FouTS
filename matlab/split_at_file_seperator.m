function fibre_objects = split_at_file_seperator(elems, file_seperator, varargin)
  

  if (size(elems, 2) ~= 3)
    error (['Dimension 2 of elements must be 3 (' num2str(size(elems,2)) ')']);
  end

  if (length(file_seperator) ~= 3)
    error (['Length of file seperator must be 3 (' length(file_seperator) ')']);
  end
  
  seperators = ones(size(elems,1),1);
  
  for i=1:length(file_seperator) 

    if isnan(file_seperator(i))
      sep_column = isnan(elems(:,i));
    elseif isinf(file_seperator(i)) && (file_seperator(i) > 0)
      sep_column = isinf(elems(:,i)) .* (elems(:,i) > 0);
    elseif isinf(file_seperator(i)) && (file_seperator(i) < 0)
      sep_column = isinf(elems(:,i)) .* (elems(:,i) < 0);
    else
      error(['Invalid file seperator, ''' num2str(file_seperator(i)) ''', must be nan, inf or -inf']);
    end
    
    seperators = seperators .* sep_column;
    
  end
  
  seperator_indices = [0; find(seperators)];

  num_fibres = size(seperator_indices,1)-1;
  
  fibre_objects = cell(num_fibres,1);
  
  for fibre_i = 1:num_fibres
    
    fibre_objects{fibre_i} = elems((seperator_indices(fibre_i)+1):(seperator_indices(fibre_i+1)-1),:);
    
  end



end