function [keys, values] = read_element_properties(location)

  disp('Reading properties...');

  file = fopen(location, 'r');
	
  if (file == -1)
    error([ 'Could not open file ' location '!' ]);
  end

     
  header_line = fgetl(file);
      
  %Seperate the header tag from the key list.
  [header, key_string] = strtok(header_line,':');
    
  if ~strcmp(header, '%%% Extended Properties File %%% - keys')
    error(['Extended properties header tag ''%%% Extended Properties File %%% - keys'' was not found at start of extended properties file ''' location '']);
  end  

  keys = parse_properties_row(key_string(3:end));
  
  num_keys = length(keys);
  
  values = cell(0, num_keys);
  
  row_i = 1;
  
  value_line = fgetl(file);
  
  while value_line ~= -1;
    
    value_row = parse_properties_row(value_line);
    
    if size(value_row,2) ~= num_keys
      error(['Number of values on line ' num2str(row_i) ' (' num2str(size(value_row,2)) ') does not match number of keys (' num2str(num_keys) ').']);
    end

    values(row_i, :) = value_row;
    
    value_line = fgetl(file);    
    row_i = row_i + 1; %#ok<NASGU>
    
  end
  
    disp('Finished reading properties...');
  
end
