function [keys, values] = read_element_properties(location)

  file = fopen(location, 'r');
	
	if (file == -1) 
		error([ 'Could not open file ' filename '!' ]);
  end

     
  header_line = fgetl(file);
      
  %Seperate the header tag from the key list.
  [header, key_string] = strtok(header_line,':');
    
  if ~strcmp(header, '%%% Extended Properties File %%% - keys')
    error(['Extended properties header tag ''%%% Extended Properties File %%% - keys'' was not found at start of extended properties file ''' location '']);
  end  

  %String ': ' from start of keys and trailing tab from end.
  key_string = deblank(key_string(3:end));

  keys = cell(0);

  keys = regexp(key_string,'\t','split')';
  
  num_keys = size(keys,1);
  
  values = cell(0, num_keys);
  
  row_i = 1;
  
  value_line = fgetl(file);
  
  while value_line ~= -1;
    
    value_row = regexp(value_line,'\t','split');
    
    %Remove the empty value after the final tab.
    value_row = value_row(1:end-1);
    
    if size(value_row,2) ~= num_keys
      error(['Number of values on line ' num2str(row_i) ' (' num2str(size(value_row,2)) ') does not match number of keys (' num2str(num_keys) ').']);
    end
    
    for key_i=1:num_keys
      values{row_i, key_i} = value_row{key_i};
    end
    
    value_line = fgetl(file);    
    row_i = row_i + 1; %#ok<NASGU>
    
  end
  
end