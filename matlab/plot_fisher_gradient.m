function plot_fisher_gradient(varargin)
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
              'coords', [],                   'matrix_:x3', 'Selected coordinates to plot.';...
              'indices', [],                  'matrix_:x:', 'Selected indices to plot.';...
              'test_include', [],             'matrix_1x:', 'Tests to include in the plots.';...                                 
             };         


  parse_arguments      
  if (help_display) 
    return;
  end

%   End arguments %


   extension = file_extension(comparison_files);
  base = file_base(comparison_files);

  if ~strcmp(extension, 'tnr')
    error(['Files must have ''tnr'' extension for Fisher gradients.']);
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
  num_tensors = size(analytic,3);
  
  if size(analytic,2) ~= num_parameters
    error('Matrices not square.');
  end
  
  if mod(num_tensors, num_parameters) ~= 0
    error(['Number of parameters (' num2str(num_parameters) ') does not divide into number of matrices (' num2str(num_tensors) ').']);
  end
  
  num_tests = num_tensors / num_parameters;
  
  if ~isempty(coords)
  
    analytic_block = reshape(analytic, [num_parameters, num_parameters, num_parameters, num_tests]);
    numeric_block = reshape(numeric, [num_parameters, num_parameters, num_parameters, num_tests]);
   
    sum_ana = sum(analytic_block,4);
    sum_num = sum(numeric_block,4);
    
    ana_bin = abs(sum_ana) > 0;
    num_bin = abs(sum_num) > 0;
    
    tog_bin = ana_bin .* num_bin;
    just_ana = ana_bin - tog_bin;
    just_num = num_bin - tog_bin;
    
    [t1,t2,t3] = ind2sub([num_parameters,num_parameters,num_parameters],find(tog_bin));
    [a1,a2,a3] = ind2sub([num_parameters,num_parameters,num_parameters],find(just_ana));
    [n1,n2,n3] = ind2sub([num_parameters,num_parameters,num_parameters],find(just_num));
    
    tog = [t1,t2,t3]
    ana = [a1,a2,a3]
    num = [n1,n2,n3]
    
    num_coords = size(coords,1);
    
    for coord_i=1:num_coords
      
      analytic = squeeze(analytic_block(coords(coord_i,1),coords(coord_i,2),coords(coord_i,3),:));
      numeric = squeeze(numeric_block(coords(coord_i,1),coords(coord_i,2),coords(coord_i,3),:));
      
      my_figure( [num2str(coords(coord_i)) ' - ' labels{coords(coord_i,1)} '::' labels{coords(coord_i,2)} '::' labels{coords(coord_i,3)} ] ,coord_i,num_coords,[]);
    
      plot(numeric, analytic, style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

      set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
      set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

      hold on;

      overall_max_value = max(max([numeric,analytic])) * 1.1;
      overall_min_value = min(min([numeric,analytic])) * 1.1;

      plot([overall_min_value; overall_max_value], [overall_min_value; overall_max_value], '-', 'color', [0.7412    0.6549    0.9608]);

      hold off;     
      
    end
    
  elseif ~isempty(indices)
    
    analytic_block = reshape(analytic, [num_parameters^3, num_tests]);
    numeric_block = reshape(numeric, [num_parameters^3, num_tests]);

    num_indices = length(indices);
    
    for index_i=1:num_indices
      
      index = indices(index_i);
      
      analytic = analytic_block(index,:);
      numeric = numeric_block(index,:);
      
      [ax1,ax2,ax3] = ind2sub([num_parameters,num_parameters,num_parameters],index);
      
      my_figure( [num2str(index) ' - ' labels{ax1} '::' labels{ax2} '::' labels{ax3} ] ,index_i,num_indices,[]);
    
      plot(numeric, analytic, style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

      set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
      set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

      hold on;

      overall_max_value = max(max([numeric,analytic])) * 1.1;
      overall_min_value = min(min([numeric,analytic])) * 1.1;

      plot([overall_min_value; overall_max_value], [overall_min_value; overall_max_value], '-', 'color', [0.7412    0.6549    0.9608]);

      hold off;     
      
    end
    
    
      
  else
        
    analytic = analytic(:);
    numeric = numeric(:);
      
    my_figure('All parameters',1,1,[]);
    
    plot(numeric, analytic, style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

    set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
    set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

    hold on;
    
    overall_max_value = max(max([numeric,analytic])) * 1.1;
    overall_min_value = min(min([numeric,analytic])) * 1.1;

    plot([overall_min_value; overall_max_value], [overall_min_value; overall_max_value], '-', 'color', [0.7412    0.6549    0.9608]);

    hold off;      
    
  end
  
  disp('')
  disp(['Plotted ' num2str(num_tests) ' tests.']);   
  disp('')
  
end