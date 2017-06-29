function plot_hessian(varargin)
%  
% PURPOSE: Plots the difference between the analytically calculated gradient and the numerically calculated gradient.
%  
% ARGUMENTS: 
%  
%           comparison_files     The filename containing the comparison between the analytical and numerical gradients.
%  
% OPTIONS (name, description, type, default):
%  
%          -style        
%                 Style of the plot.
%                 string
%                 'x'
%  
%          -marker_size  
%                 Size of marker
%                 int
%                 10
%  
%          -line_width   
%                 Width of line
%                 int
%                 2
%  
%          -param_include
%                 Parameters to include in the plots.
%                 matrix_1x:
%  
%          -test_include
%                 Tests to include in the plots.
%                 matrix_1x:
                
                
  description = 'Plots the difference between the analytically calculated gradient and the numerically calculated gradient.';
  
  arguments = {'comparison_files', 'The filename containing the comparison between the analytical and numerical gradients.'};
  
  options = {...
              'style        ',  'x',           'string',     'Style of the plot.';...                       
              'marker_size  ',  10,            'int',        'Size of marker';...  
              'line_width   ',  2,             'int',        'Width of line';...                
              'index_include',  [],            'matrix_1x:', 'Parameter indices to include in the plots (overrides coord_include).';...                                
              'coord_include',  [],            'matrix_:x:', 'Parameter coordinates to include in the plots.';...  
              'test_include',   [],            'matrix_1x:', 'Tests to include in the plots.';...       
              'diagonal',       0,             'bool',       'Only displays diagonal components';... 
              'collapse',       0,             'bool',       'Collapse all plots into single one';...
            };         


  num_args = nargin;
  parse_arguments      
  if (help_display) 
    return;
  end

%   End arguments %


  extension = file_extension(comparison_files);
  base = file_base(comparison_files);

  if ~strcmp(extension, 'tnr')
    error(['Files must have ''tnr'' extension for tensors.']);
  end
  
  base_extension = file_extension(base);
  base_base = file_base(base);

  [analytic, labels] = load_tensor([base_base '.analytic.' base_extension '.' extension]);
  numeric = load_tensor([base_base '.numeric.' base_extension '.' extension]);    

  if any(size(analytic) ~= size(numeric))
    error(['Size of analytic (' mat2str(size(analytic)) ') does not match size of numeric (' mat2str(size(numeric)) ').']);
  end
  
  %Ensure matrices are of the same size.
  num_parameters = size(analytic,1);
  num_tests = size(analytic,3);
  
  if size(analytic,2) ~= num_parameters
    error('Matrices not square.');
  end
  
  %Set to display to all comparisons if not selection explicitly provided.
  if diagonal
    for param_i = 1:num_parameters
      coord_include(param_i,:) = [param_i param_i];
    end
    
  elseif ~isempty(index_include)
    
    [row col] = ind2sub([num_parameters, num_parameters], index_include);
    
    coord_include = [row' col'];
    
  elseif isempty(coord_include)
    param_count = 1;
    for row_i = 1:num_parameters
      for col_i = 1:num_parameters
        coord_include(param_count,:) = [row_i col_i];
        param_count = param_count +1;
      end
    end
  end
  
  if isempty(test_include)
    test_include = 1:num_tests;
  end
  
  num_param_include = size(coord_include,1);
  
  if collapse
    my_figure('All parameters');
  end
  
  %Plot comparisons.
  for param_index_i = 1:num_param_include
    
    param_index = coord_include(param_index_i,:);
      
    numeric_param = squeeze(numeric(param_index(1), param_index(2), test_include));
    analytic_param = squeeze(analytic(param_index(1), param_index(2), test_include));
    
    if ~collapse
      my_figure([labels{param_index(1)}(1:end) ' :: ' labels{param_index(2)}(1:end)], param_index_i, num_param_include);
    end
    
    plot(numeric_param, analytic_param, style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

    set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
    set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

    min_value = min([numeric_param; analytic_param]) * 1.1;
    max_value = max([numeric_param; analytic_param]) * 1.1;

    hold on;

    plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

    hold off;      

  end  

  disp('')
  disp(['Plotted ' num2str(num_tests) ' tests.']);   
  disp('')
  
end
