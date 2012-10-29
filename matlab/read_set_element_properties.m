function [keys, values] = read_set_element_properties(location)

  file = fopen(location, 'r');
	
	if (file == -1) 
    
    keys = cell(0);
    values = cell(0);
    return;
  end

     
  header_line = fgetl(file);
      
  %Seperate the header tag from the key list.
  [header, key_string] = strtok(header_line,':');
    
  if ~strcmp(header, '%%% Extended Properties File %%% - keys')
    error(['Extended properties header tag ''%%% Extended Properties File %%% - keys'' was not found at start of extended properties file ''' location '']);
  end  

  %String ': ' from start of keys and trailing tabs from end.
  keys = parse_properties_row(key_string(3:end));
  
  num_keys = length(keys);

  row_i = 1;
  
  value_line = fgetl(file);
  
  while value_line ~= -1;
       
    sub_row_i = 1;
    
    while length(value_line) < 8 || ~strcmp(value_line(1:8), '--- END ') || ~strcmp(value_line((end-3):end), ' ---')  
    
      value_row = parse_properties_row(value_line);

      if size(value_row,2) ~= num_keys
        error(['Number of values on line ' num2str(row_i) ' (' num2str(size(value_row,2)) ') does not match number of keys (' num2str(num_keys) ').']);
      end

      values{row_i,1}(sub_row_i,:) = value_row;


      value_line = fgetl(file);    
      sub_row_i = sub_row_i + 1; %#ok<NASGU>
      
      if value_line == -1
        error (['Found end of file before row delimeter for row (in row ' num2str(row_i) ').']);
      end      
      
    end

    row_i = row_i + 1; %#ok<NASGU>
    value_line = fgetl(file);
        
  end
  
  
  
end