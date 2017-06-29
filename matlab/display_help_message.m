function display_help_message(description, arguments, options)

disp(' ');

disp(['PURPOSE: ' description]);

disp(' ');

disp(['ARGUMENTS: ']);

for arg_i = 1:size(arguments,1)
  
  disp(' ');
  disp(['          ' arguments{arg_i,1} '     ' arguments{arg_i,2}]);
  
end  

disp(' ');

disp('OPTIONS (name, description, type, default):');

for option_i = 1:size(options,1)
  
  disp(' ');
  disp(['         -' options{option_i,1}]);
  disp(['                ' options{option_i,4}]);
  disp(['                ' options{option_i,3}]);
  
  if strcmp(options{option_i,3}, 'string')
    disp(['                ''' options{option_i,2} '''' ]);    
  else

    if ~isstr(options{option_i,2})
      options{option_i,2} = mat2str(options{option_i,2});
    end
    
    if size(strfind(options{option_i,3}, 'matrix'),1) == 0
      disp(['                ' options{option_i,2} ]);    
    end
  end

end

disp(' ');
