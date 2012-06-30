function plot_gradient(varargin)
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
              'style        ', 'x',           'string',     'Style of the plot.';...                       
              'marker_size  ', 10,            'int',        'Size of marker';...  
              'line_width   ', 2,             'int',        'Width of line';...                
              'param_include', [],            'matrix_1x:', 'Parameters to include in the plots.';...                                
              'test_include', [],             'matrix_1x:', 'Tests to include in the plots.';...                                 
            };         


  parse_arguments      
  if (help_display) 
    return;
  end

%   End arguments %


  extension = file_extension(comparison_files);
  base = file_base(comparison_files);

  %Load numeric and analytic gradients in 'unzipped' (flattened into vectors) formats.
  if strcmp(extension, 'str')
    
    [analytic, labels] = load_unzip_strands([base '.analytic.' extension]);
    numeric = load_unzip_strands([base '.numeric.' extension]);
  
  elseif strcmp(extension, 'tct')
    
    [analytic, labels] = load_unzip_tracts([base '.analytic.tct']);
    numeric = load_unzip_tracts([base '.numeric.' extension]);
  
  elseif strcmp(extension, 'sst')

    [analytic, labels] = load_unzip_strand_sets([base '.analytic.' extension]);
    numeric = load_unzip_strand_sets([base '.numeric.' extension]);

  elseif strcmp(extension, 'tst')

    [analytic, labels] = load_unzip_tract_sets([base '.analytic.' extension]);
    numeric = load_unzip_tract_sets([base '.numeric.' extension]);

  elseif strcmp(extension, 'ssc')

    [analytic, labels] = load_unzip_strand_sections([base '.analytic.' extension]);
    numeric = load_unzip_strand_sections([base '.numeric.' extension]);

  elseif strcmp(extension, 'tsc')

    [analytic, labels] = load_unzip_tract_sections([base '.analytic.' extension]);
    numeric = load_unzip_tract_sections([base '.numeric.' extension]);

  elseif strcmp(extension, 'trp')

    [analytic, labels] = load_unzip_triples([base '.analytic.' extension]);
    numeric = load_unzip_triples([base '.numeric.' extension]);
    
  elseif strcmp(extension, 'mif')

    [analytic, labels] = load_unzip_image([base '.analytic.' extension]);
    numeric = load_unzip_image([base '.numeric.' extension]);

  elseif strcmp(extension, 'sta')
    
    [analytic, labels] = load_unzip_mcmc_state([base '.analytic.' extension]);
    numeric = load_unzip_mcmc_state([base '.numeric.' extension]);
    
  else
    error(['Unrecognised extension ''' char(extension) '''.']);
  end 
  
  
  %Ensure matrices are of the same size.
  num_parameters = size(analytic,1);
    
  if num_parameters ~= size(numeric,1);
    error (['Number of parameters in numeric matrix (' num2str(size(numeric,1)) ') does not match number in analytic matrix (' num2str(num_parameters) ').']);
  end
  
  num_tests = size(analytic,2);
  
  if num_tests ~= size(numeric,2);
    error (['Number of parameters in numeric matrix (' num2str(size(numeric,2)) ') does not match number in analytic matrix (' num2str(num_tests) ').']);
  end
  
  %Set to display to all comparisons if not selection explicitly provided.
  if isempty(param_include)
    param_include = 1:num_parameters;
  end
  
  if isempty(test_include)
    test_include = 1:num_tests;
  end
  
  num_param_include = length(param_include(:));
  
  %Plot comparisons.
  for param_index_i = 1:num_param_include
    
    param_index = param_include(param_index_i);
      
    numeric_param = numeric(param_index, test_include)';
    analytic_param = analytic(param_index, test_include)';
    
    my_figure(labels{param_index}, param_index_i, num_param_include);
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
