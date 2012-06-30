function [keys, values] = load_extend_elem_props(location)

  keys = cell(0);
  values = [];

  if exist(location)
    
    ext_file = fopen(location, 'r');
    header_line = fgetl(ext_file);
    fclose(ext_file);
    
    %Seperate the header tag from the key list.
    [header, key_string] = strtok(header_line,':');
    
    if ~strcmp(header, '%%% Extended Properties File %%% - keys')
      error(['Extended properties header tag ''%%% Extended Properties File %%% - keys'' was not found at start of extended properties file ''' location '']);
    end  
     
    %String ': ' from start of keys and trailing tab from end.
    key_string = deblank(key_string(3:end));
    keys = cell(0);
    
    while ~isempty(key_string)
      [keys{end+1}, key_string] = strtok(key_string);
    
    end
    
    values = load(location);
    
  end
  
end