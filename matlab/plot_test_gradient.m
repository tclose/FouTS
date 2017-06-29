function plot_test_gradient(varargin)
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
%          -include      
%                 Samples to include in the plots.
%                 matrix_1x:
                
                
  description = 'Plots the difference between the analytically calculated gradient and the numerically calculated gradient.';
  
  arguments = {'comparison_files', 'The filename containing the comparison between the analytical and numerical gradients.'};
  
  options = {...
              'style        ', 'x',           'string',     'Style of the plot.';...                       
              'marker_size  ', 10,            'int',        'Size of marker';...  
              'line_width   ', 2,             'int',        'Width of line';...                
              'include      ', [],            'matrix_1x:', 'Samples to include in the plots.';...                                
            };         


  num_args = nargin;
  parse_arguments      
  if (help_display) 
    return;
  end

%   End arguments %


  file_parts = regexp(comparison_files, '\.', 'split');

  extension = file_parts(end);
  base = cat(2,file_parts{1:end-1});
  

  
  if strcmp(extension, 'str')
    
    analytic = load_strands([base '.analytic.str']);
    numeric = load_strands([base '.numeric.str']);

    num_strands = size(analytic,1);
    
    if num_strands == 0
      error ('No strands loaded.');
    end
    
    num_degree = length(analytic{1});
    
    if num_strands ~= size(numeric,1)
      error(['Size of numeric gradient (' num2str(size(numeric,1)) ') did not match size of analytic gradient (' num2str(num_strands) ').']);
    end

    if isempty(include)
      include = 1:num_strands;
    end    
    
    num_include = length(include);
    
    degrees = cell(num_degree);
    
    for degree_i = 1:num_degree
      degrees{degree_i} = zeros(3*num_include,2);
    end
    
    for strand_i = include
      
      analytic_strand = analytic{strand_i};
      numeric_strand = numeric{strand_i};
      
      start_index = (strand_i-1)*3 + 1;
      end_index = start_index+2;
      
      for degree_i = 1:num_degree

        degrees{degree_i}(start_index:end_index,1) = numeric_strand(degree_i,:)';
        degrees{degree_i}(start_index:end_index,2) = analytic_strand(degree_i,:)';        

        
      end  
      
      
    end
    

    for degree_i = 1:num_degree
        
      degree = degrees{degree_i};
      
      my_figure([comparison_files '- degree ' num2str(degree_i-1)], degree_i, num_degree, []);
      h = plot(degree(:,1), degree(:,2), style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

      set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
      set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

      min_value = min(min(degree)) * 1.1;
      max_value = max(max(degree)) * 1.1;
      
      hold on;

      plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

      hold off;      
      
    end  
      
    
  elseif strcmp(extension, 'sst')
    
    analytic = load_strand_sets([base '.analytic.sst']);
    numeric = load_strand_sets([base '.numeric.sst']);

    num_sets = size(analytic,1);
    
    if num_sets == 0
      error('No sets loaded.');
    elseif size(analytic{1},1) == 0
      error ('Loaded have size 0.');
    end
    
    num_degree = size(analytic{1}{1},1);
    
    if num_sets ~= size(numeric,1)
      error(['Size of numeric gradient (' num2str(size(numeric,1)) ') did not match size of analytic gradient (' num2str(num_sets) ').']);
    end

    if isempty(include)
      include = 1:num_sets;
    end    
    
    num_include = length(include);
    
    degrees = cell(num_degree,1);
    
    for degree_i = 1:num_degree
      degrees{degree_i} = zeros(3*num_include,2);
    end
    
    
    start_index = 1;
    
    for set_i = include
      
      analytic_set = analytic{set_i};
      numeric_set = numeric{set_i}; 
        
      num_strands = size(analytic_set,1);
      
      if (num_strands ~= size(numeric_set,1))
        error (['Size of numeric gradient, ' num2str(size(numeric_set,1)) ' in set ' num2str(set_i) ' does not match analytic gradient, ' num2str(num_strands) '.']);
      end
      
      for strand_i = 1:num_strands

          analytic_strand = analytic_set{strand_i};
          numeric_strand = numeric_set{strand_i};

          for degree_i = 1:num_degree

            degrees{degree_i}(start_index:(start_index+2),1) = numeric_strand(degree_i,:)';
            degrees{degree_i}(start_index:(start_index+2),2) = analytic_strand(degree_i,:)';        

          end  
          
          start_index = start_index + 3;
      end
      
    end
    
    size(degrees{1})

    for degree_i = 1:num_degree
        
      degree = degrees{degree_i};
      
      my_figure([comparison_files '- degree ' num2str(degree_i-1)], degree_i, num_degree, []);
      h = plot(degree(:,1), degree(:,2), style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

      set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
      set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

      min_value = min(min(degree)) * 1.1;
      max_value = max(max(degree)) * 1.1;
      
      hold on;

      plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

      hold off;      
      
    end  
      
  elseif strcmp(extension, 'tst')
    
    analytic = load_tract_sets([base '.analytic.tst']);
    numeric = load_tract_sets([base '.numeric.tst']);

    num_sets = size(analytic,1);
    
    if num_sets == 0
      error('No sets loaded.');
    elseif size(analytic{1},1) == 0
      error ('Loaded have size 0.');
    end
    
    num_degree = size(analytic{1}{1},1);
    
    if num_sets ~= size(numeric,1)
      error(['Size of numeric gradient (' num2str(size(numeric,1)) ') did not match size of analytic gradient (' num2str(num_sets) ').']);
    end

    if isempty(include)
      include = 1:num_sets;
    end    
    
    num_include = length(include);
    
    degrees = cell([3,num_degree]);
    
    for strand_i = 1:3
      for degree_i = 1:num_degree
        degrees{strand_i, degree_i} = zeros(3*num_include,2);
      end
    end
    
    start_index = 1;
    
    for set_i = include
      
      analytic_set = analytic{set_i};
      numeric_set = numeric{set_i}; 
        
%      num_strands = size(analytic_set,1);
      
      
      for strand_i = 1:3

          analytic_strand = analytic_set{strand_i};
          numeric_strand = numeric_set{strand_i};

          for degree_i = 1:num_degree

            degrees{strand_i,degree_i}(start_index:(start_index+2),1) = numeric_strand(degree_i,:)';
            degrees{strand_i,degree_i}(start_index:(start_index+2),2) = analytic_strand(degree_i,:)';        

          end  
          
          start_index = start_index + 3;
      end
      
    end
    
    size(degrees{1})

    for strand_i = 1:3
      for degree_i = 1:num_degree

        degree = degrees{strand_i, degree_i};

        my_figure([comparison_files '- axis ' num2str(strand_i) ' degree ' num2str(degree_i-1)], (strand_i-1)*num_degree + degree_i, 3 * num_degree, []);
        h = plot(degree(:,1), degree(:,2), style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

        set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
        set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

        min_value = min(min(degree)) * 1.1;
        max_value = max(max(degree)) * 1.1;

        hold on;

        plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

        hold off;      

      end  
    end
    
  elseif strcmp(extension, 'tct')

    analytic = load_tracts([base '.analytic.tct']);
    numeric = load_tracts([base '.numeric.tct']);

    num_tracts = size(analytic,1);
    
    if num_tracts == 0
      error ('No tracts loaded.');
    end
    
    num_degree = length(analytic{1,1});
    
    if num_tracts ~= size(numeric,1)
      error(['Size of numeric gradient (' num2str(size(numeric,1)) ') did not match size of analytic gradient (' num2str(num_tracts) ').']);
    end

    if isempty(include)
      include = 1:num_tracts;
    end    
    
    num_include = length(include);
    
    axis_degrees = cell(3,num_degree);
    
    for axis_i = 1:3
      for degree_i = 1:num_degree
        axis_degrees{axis_i, degree_i} = zeros(3*num_include,2);
      end
    end
    
    for tract_i = include

      start_index = (tract_i-1)*3 + 1;
      end_index = start_index+2;      
      
      for axis_i = 1:3
        
        analytic_strand = analytic{tract_i, axis_i};
        numeric_strand = numeric{tract_i, axis_i};


        for degree_i = 1:num_degree

          axis_degrees{axis_i, degree_i}(start_index:end_index,1) = analytic_strand(degree_i,:)';
          axis_degrees{axis_i, degree_i}(start_index:end_index,2) = numeric_strand(degree_i,:)';        


        end  
        
      end
      
    end
    
    for axis_i = 1:3
      for degree_i = 1:num_degree

        degree = axis_degrees{axis_i, degree_i};

        my_figure([comparison_files '- Axis ' num2str(axis_i) ', Degree ' num2str(degree_i)], (axis_i-1)*num_degree + degree_i, 3 * num_degree, []);
        h = plot(degree(:,1), degree(:,2), style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

        set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
        set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

        min_value = min(min(degree)) * 1.1;
        max_value = max(max(degree)) * 1.1;

        hold on;

        plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

        hold off;      

      end  
    end
    
    
  elseif strcmp(extension, 'trp')
    
    analytic = load_triples([base '.analytic.trp']);
    numeric = load_triples([base '.numeric.trp']);
    

    num_triples = size(analytic,1);
    
    if num_triples ~= size(numeric,1)
      error(['Size of numeric gradient (' num2str(size(numeric,1)) ') did not match size of analytic gradient (' num2str(num_triples) ').']);
    end    
    
    
    if isempty(include)
      include = 1:num_triples;
    end    
    
    num_include = length(include);
    
    for dim_i = 1:3
      
      my_figure([comparison_files ' - Dim ' num2str(dim_i) '.'], dim_i, 3, []);
      
      h = plot(analytic(include,dim_i), numeric(include,dim_i), style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

      set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
      set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

      min_value = min(min([ analytic(include,dim_i) numeric(include,dim_i)])) * 1.1;
      max_value = max(max([ analytic(include,dim_i) numeric(include,dim_i)])) * 1.1;

      hold on;

      plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

      hold off;    
      
    end

  elseif strcmp(extension, 'ssc')    
    
    analytic = load_strand_sections([base '.analytic.ssc']);
    numeric = load_strand_sections([base '.numeric.ssc']);
    

    num_sections = size(analytic,1);
    
    if num_sections ~= size(numeric,1)
      error(['Size of numeric gradient (' num2str(size(numeric,1)) ') did not match size of analytic gradient (' num2str(num_segments) ').']);
    end    
    
        
    if isempty(include)
      include = 1:num_sections;
    end    
    
    num_include = length(include);
    
    
    position_and_tangent = cell(2);
    
    for pos_or_tang = 1:2
      position_and_tangent{pos_or_tang} = zeros(3*num_include,2);
    end
    
    for section_i = include
      
      start_index = (section_i-1)*3 + 1;
      end_index = start_index+2;
      
      for pos_or_tang = 1:2

        position_and_tangent{pos_or_tang}(start_index:end_index,1) = numeric{section_i}(pos_or_tang,:)';
        position_and_tangent{pos_or_tang}(start_index:end_index,2) = analytic{section_i}(pos_or_tang,:)';        

        
      end  
      
      
    end
    

    
    for pos_or_tang = 1:2
      
      position_or_tangent = position_and_tangent{pos_or_tang};
      
      if (pos_or_tang == 1)
        title_string = [comparison_files ' - Position.'];
      else
        title_string = [comparison_files ' - Tangent.'];
      end
      
      my_figure(title_string, pos_or_tang, 2, []);
      
      h = plot(position_or_tangent(:,1), position_or_tangent(:,2), style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

      set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
      set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

      min_value = min(min(position_or_tangent)) * 1.1;
      max_value = max(max(position_or_tangent)) * 1.1;

      hold on;

      plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

      hold off;    
      
    end
  
  elseif strcmp(extension, 'tsc')    
    
    analytic = load_tract_sections([base '.analytic.tsc']);
    numeric = load_tract_sections([base '.numeric.tsc']);
    

    num_sections = size(analytic,1);
    
    if num_sections ~= size(numeric,1)
      error(['Size of numeric gradient (' num2str(size(numeric,1)) ') did not match size of analytic gradient (' num2str(num_sections) ').']);
    end    
    
        
    if isempty(include)
      include = 1:num_sections;
    end    
    
    num_include = length(include);
    
    
    section_components = cell(4);
    
    for comp_i = 1:4
      section_components{comp_i} = zeros(3*num_include,2);
    end
    
    for section_i = include
      
      start_index = (section_i-1)*3 + 1;
      end_index = start_index+2;
      
      for comp_i = 1:4

        section_components{comp_i}(start_index:end_index,1) = numeric{section_i}(comp_i,:)';
        section_components{comp_i}(start_index:end_index,2) = analytic{section_i}(comp_i,:)';        

      end  
      
      
    end
    

    
    for comp_i = 1:4
      
      component = section_components{comp_i};
      
      if (comp_i == 1)
        title_string = [comparison_files ' - Position.'];
      elseif comp_i == 2
        title_string = [comparison_files ' - Tangent.'];
      elseif comp_i == 3
        title_string = [comparison_files ' - Width1.'];
      else
        title_string = [comparison_files ' - Width2.'];
      end
      
      my_figure(title_string, comp_i, 4, []);
      
      h = plot(component(:,1), component(:,2), style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

      set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
      set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

      min_value = min(min(component)) * 1.1;
      max_value = max(max(component)) * 1.1;

      hold on;

      plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

      hold off;    
      
    end
          
    
  elseif strcmp(extension, 'mif')
    
    analytic_struct = read_image([base '.analytic.mif']);
    numeric_struct  = read_image([base '.numeric.mif']);
    
    analytic = analytic_struct.data;
    numeric = numeric_struct.data;
    
    num_include = size(analytic,1) * size(analytic,2) * size(analytic,3) * size(analytic,4);
    
    analytic = reshape(analytic, num_include, 1);
    numeric = reshape(numeric, num_include, 1);
    
    min_value = min(min([ analytic numeric])) * 1.1;
    max_value = max(max([ analytic numeric])) * 1.1;
    
    my_figure([comparison_files], []);

    h = plot(analytic, numeric, style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

    set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
    set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

    min_value = min(min([ analytic numeric])) * 1.1;
    max_value = max(max([ analytic numeric])) * 1.1;

    hold on;

    plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

    hold off;  

  elseif strcmp(extension, 'tes')
    
    analytic = load([base '.analytic.tes']);
    numeric  = load([base '.numeric.tes']);
    
    num_include = size(analytic,1);
    
    
    for dim_i = 1:2
    
      min_value = min(min([analytic(:,dim_i) numeric(:,dim_i)])) * 1.1;
      max_value = max(max([analytic(:,dim_i) numeric(:,dim_i)])) * 1.1;

      my_figure([comparison_files ' - dim ' num2str(dim_i)], dim_i, 2, []);

      h = plot(analytic(:,dim_i), numeric(:,dim_i), style, 'MarkerSize', marker_size, 'LineWidth',line_width);    

      set(get(gca,'xlabel'), 'string', 'numeric', 'color', [1 1 1]);
      set(get(gca,'ylabel'), 'string', 'analytic', 'color', [1 1 1]);   

      min_value = min(min([ analytic(:,dim_i) numeric(:,dim_i)])) * 1.1;
      max_value = max(max([ analytic(:,dim_i) numeric(:,dim_i)])) * 1.1;

      hold on;

      plot([min_value; max_value], [min_value; max_value], '-', 'color', [0.7412    0.6549    0.9608]);

      hold off;    
    end
    
  else
    error(['Unrecognised extension ''' char(extension) '''.']);
  end 


  fprintf('\n');
  disp(['Plotted ' num2str(num_include) ' comparisons.']);   
  
  
end
