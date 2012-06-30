function [tensor_matrix, labels] = load_tensor(filename)
  
  file = fopen(filename, 'r');
	
  if (file == -1)
    error([ 'Could not open file ' filename '!' ]);
  end
     
  header_line = fgetl(file);
  
  fclose(file);
  
    %Seperate the header tag from the key list.
  [preamble, header] = strtok(header_line,':');
    
  if ~strcmp(preamble, '%%% Keys %%%')
    error(['Did not find tensor preamble at location ''' filename '']);
  end  
  
  header = deblank(header(2:end));
  
  if isempty(header)
    error('No keys found in tensor file.');
  elseif isempty(regexp(header,' ', 'once'))
    labels{1} = line;
  else  
    labels = regexp(header,' ','split');
  end
  
  tensor_matrix = load(filename)';
  
  num_elements = size(tensor_matrix,1);
  
  num_objects = size(tensor_matrix,2)/num_elements;
  
   
  if num_objects == 0
    error('No elements loaded from file.')
  elseif floor(num_objects) ~= num_objects
    error(['Number of elements (' num2str(num_elements) ') does not divide into number of rows (' num2str(size(tensor_matrix,2)) ').']);
  end
    
  if num_elements ~= length(labels)
    error(['Number of elements (' num2str(num_elements) ') does not match number of labels (' num2str(length(labels)) ').']);
  end
  
  tensor_matrix = reshape(tensor_matrix, num_elements, num_elements, num_objects);
  
end