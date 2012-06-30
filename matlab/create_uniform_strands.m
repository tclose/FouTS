function create_uniform_strands(output_dir, varargin)
%  
% PURPOSE: Initialise a number of straight strands evenly distributed around a sphere.
%  
% ARGUMENTS: 
%  
%           output_directory     The directory that the newly created strands will be saved to
%  
% OPTIONS (name, description, type, default):
%  
%          -num_strands   
%                 number of strands to be created
%                 int
%                 16
%  
%          -strand_length 
%                 The v1 magnitude of the fourier descriptors of the strands
%                 float
%                 0.053         
                

  description = 'Initialise a number of straight strands evenly distributed around a sphere.';
  
  arguments = {'output_directory', 'The directory that the newly created strands will be saved to'};

  options = {...
            'num_strands   ', 16,        'int',    'number of strands to be created';...
            'strand_length ', 0.2755,    'float',  'The v1 magnitude of the fourier descriptors of the strands';...            
            };         
                                    
  supplied_options = parse_options(options, varargin);

  for option_i = 1:size(supplied_options,1)
    if strcmp(supplied_options{option_i,3}, 'string')
      eval([supplied_options{option_i,1} ' = ''' supplied_options{option_i,2} ''';']);
    else
      
      if ~isstr(supplied_options{option_i,2})
        supplied_options{option_i,2} = mat2str(supplied_options{option_i,2});
      end
      
      eval([supplied_options{option_i,1} ' = ' supplied_options{option_i,2} ';']);
    end
  end


  if help_display
    display_help_message(description, arguments, options);
    return;
  end
  
  
  degree_1s = strand_length * equidistribute(num_strands*2);
  
  degree_1s_z = degree_1s(:,3);
  
  positive_degree_1s_z = find(degree_1s_z >= 0);
  
  degree_1s = degree_1s(positive_degree_1s_z,:);
  
  
  for strand_i = 1:num_strands
    
    strand = [0 0 0; degree_1s(strand_i,:); 0 0 0];
    
    
    save([output_dir '/strand_' num2str(strand_i) '.frr.txt'], 'strand', '-ASCII');
    
  end

  
  
end