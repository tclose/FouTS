function options = parse_options(options, supplied_options)
% function options = parse_options(options, supplied_options)

%   if (size(options,2) ~= 3)
%     error(['Options supplied by source file have not been formatted properly they require three values for each line: key, default_value, type_string.']);
%   end
% 
%   for option_i = 1:size(options,1)
%     
%     key = options{option_i,1};
%     default_value = options{option_i,2};
%     type = options{option_i,3};
%     descript = options{option_i,4}
%     
%     if ~isstr(key)
% 			error(['In source file, option key ' num2str(option_i) ': ' num2str(key) ' needs to be a string']);
%     end
%     
%     if ~isstr(type)
% 			error(['In source file, option type ' num2str(option_i) ': ' num2str(key) ' needs to be a string']);
%     end
%     
%     if ~isstr(descript)
% 			error(['In source file, option type ' num2str(option_i) ': ' num2str(key) ' needs to be a string']);
%     end
%     
%     check_type(default_value, type);
%     
%   end
  

  options{end+1,1} = 'help_display';
  options{end,2} = 0;
  options{end,3} = 'bool';
  options{end,4} = 'Displays the synopsis of the matlab function.';
  
 
  supplied_arg_i = 1;
  
  while supplied_arg_i <= length(supplied_options)
    
    
    option = supplied_options{supplied_arg_i};
    
    if ~isempty(option)
      key = strtrim(option);
    else
      key = [];
    end

    if ~is_key(key)
			error(['Option ' key ' is not a valid key (NB: needs to begin with the character ''-'').']); 
    end
    
    
	key = strtrim(key(2:end));
		
    key_index = [];
    
    for option_i = 1:size(options,1)
  		if strfind(options{option_i,1}, key) == 1
        key_index = [key_index; option_i];
      end
    end

	if length(key_index) == 0
		error(['Option ' key ' not a valid option.']);
	end

	if length(key_index) > 1
		error(['Option ' key ' is ambiguous.']);
    end

    value = 1;
    
    if supplied_arg_i < length(supplied_options)
      
       if ~is_key(supplied_options{supplied_arg_i+1})% Not sure why this is here && ~isempty(supplied_options{supplied_arg_i+1})
          
       		value = supplied_options{supplied_arg_i+1};
          
          % \ is used as an escape character if the value string begins with '-'.
          if ~isempty(value)
            if value(1) == '\'
              value = value(2:end);
            end
          end
          
          supplied_arg_i = supplied_arg_i + 1;
       end
    end

		type = options{key_index,3};
    
		check_type(key, value, type);

		options{key_index,2} = value;
    
    supplied_arg_i = supplied_arg_i + 1;

  end						
	
  
end			 			  