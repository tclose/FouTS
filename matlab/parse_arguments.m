function args = parse_arguments(description, arguments, options, varargin)
    % A script to automatically parse input arguments and options, and display help messages when given the '-help' option
    % requires that the following cell variables have been set in the given format and the function accepts 'varargin' as its only
    % input.
    % 
    %   description = 'Description of what the function does goes here';
    %   
    %   arguments = {'argument_name_in_function', 'Description of the argument goes here'};
    % 
    %   options = {...
    %             'argument_name_in_function  ', default_value,        'option type',  'Description of the argument goes here';...
    %                                           ... repeat for other options ...
    %             };  
    %
    % NB: Note that different from a typical bash command, the arguments must
    % appear at the start of the argument list, i.e. before the options.


    if ~exist('multiple_last_arguments','var')
      multiple_last_arguments = 0;
    end

    % Check for first argument and help option. 
    num_supplied_arguments = 0;
    
    num_args = nargin - 3;

    for arg_i = 1:num_args

      if size(varargin{arg_i},1) == 1 
        if strfind('-help_display', varargin{arg_i}) == 1
          display_help_message(description, arguments, options);
          args = 'help';
          return;
        end


        %Find the index of the first option.  This signifies the end of the
        %arguments.

        dash_indices = strfind(varargin{arg_i}, '-');

        if ~isempty(dash_indices) 
          if dash_indices(1) == 1 && ~num_supplied_arguments
            num_supplied_arguments = arg_i - 1;
          end
        end

      end
    end

    if ~num_supplied_arguments
      num_supplied_arguments = num_args;
    end

    num_arguments = size(arguments,1);

    if num_supplied_arguments < num_arguments
      error (['Not enough arguments were supplied, (found ' num2str(num_supplied_arguments) ', expected ' num2str(num_arguments) '). \n NB: Note that all arguments must be supplied before the first option.']);
    end

    if ~multiple_last_arguments && (num_supplied_arguments > num_arguments)
      error (['Too many arguments were supplied, (found ' num2str(num_supplied_arguments) ', expected ' num2str(num_arguments) ').']);
    end

    % Parse Arguments.

    if multiple_last_arguments

      for arg_i = 1:(num_arguments-1)
        eval([ arguments{arg_i,1} ' = ' mat2str(varargin{arg_i}) ';' ]);
      end

      eval([ arguments{num_arguments,1} '=' 'cell(0);']);

      for arg_i = num_arguments:num_supplied_arguments
        eval([ arguments{num_arguments,1} '{' num2str(arg_i - num_arguments + 1) '} = ' mat2str(varargin{arg_i}) ';' ]);
      end

    else

      for arg_i = 1:num_arguments
        eval([ arguments{arg_i,1} ' = ' mat2str(varargin{arg_i}) ';' ]);
      end

    end  

    % Parse Options.
    option_args = varargin(num_supplied_arguments+1:end);

    supplied_options = parse_options(options, option_args);

    for option_i = 1:size(supplied_options,1)
      if strcmp(supplied_options{option_i,3}, 'string')
        eval([supplied_options{option_i,1} ' = ''' supplied_options{option_i,2} ''';']);
      else

        if ~isstr(supplied_options{option_i,2})
          if isempty(supplied_options{option_i,2})
            supplied_options{option_i,2} = '[]';
          else
            supplied_options{option_i,2} = mat2str(supplied_options{option_i,2});
          end
        end

        eval([supplied_options{option_i,1} ' = ' supplied_options{option_i,2} ';']);
      end
    end
end

