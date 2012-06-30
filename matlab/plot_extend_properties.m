function fig = plot_extend_elem_properties(figure_name, prop_keys, prop_values, include, fig_index, num_figures)


  if ~exist('include', 'var') 
    include = [];
  end
    
  if ~exist('fig_index', 'var') 
    fig_index = 2;
  end
   
  if ~exist('num_figures', 'var') 
    num_figures = 3;
  end
      

  if ~isempty(prop_keys)
    num_values = size(prop_values,1);
     
    if isempty(include)
      include = [1:num_values];
    end

    fig = my_figure(figure_name,fig_index,num_figures,[]);

    all_props = [];

    legend_keys = cell(0);

    for key = prop_keys
      
%       if ~strcmp(key,'elapsed_time')

        props = get_properties(prop_keys, prop_values, key);
        props = props(include);

        max_abs_value = max([abs(max(props)) abs(min(props))]);

        if (max_abs_value > 0)
          sci_notation = ceil(log10(max_abs_value));
        else
          sci_notation = 1;
        end

        props = props ./ (10^sci_notation);

        all_props = [all_props, props];
        legend_keys(end+1) = strcat(strrep(key,'_',' '), ' x 10^{', num2str(sci_notation),'}');
      
%       end
      
    end
    
    set(gca, 'NextPlot', 'ReplaceChildren');
    set(gca, 'LineStyleOrder', {'-',':','--'});
    plot(include-1, all_props);
    set(gca, 'ytick', [-1:0.05:1]);
    legend(legend_keys);
  else
    disp('No extended properties found.');
    fig = [];
  end
end
