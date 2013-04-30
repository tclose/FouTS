function fig = plot_extend_elem_properties(figure_name, prop_keys, prop_values, include, fig_index, num_figures, offset_values)


  if ~exist('include', 'var') 
    include = [];
  end
    
  if ~exist('fig_index', 'var') 
    fig_index = 2;
  end
   
  if ~exist('num_figures', 'var') 
    num_figures = 3;
  end
      
  if ~exist('offset_values', 'var')
    offset_values = 1;
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

        props = get_properties(prop_keys, prop_values, key);
        props = props(include);

        min_value = min(props);
        max_value = max(props);    
        
        if ~offset_values || (max_value > 0) && (min_value < 0)
            max_abs_value = max([abs(max_value) abs(min_value)]);
            sci_not = sci_notation(max_abs_value);
            offset = 0;
            props = props ./ (10^sci_not);
        else
            sci_not = sci_notation(max_value - min_value);
            offset_exp = 10 ^ (sci_not + 1);
            
            % Check if the values don't cross zero
            if min_value > 0
                offset = floor(min_value / offset_exp) * offset_exp;
            else
                offset = ceil(max_value / offset_exp) * offset_exp;  
            end
            props = (props - offset) ./ (10^sci_not);
        end
        all_props = [all_props, props];
        legend_key = strcat(strrep(key,'_',' '), ...
            ' x 10^{', num2str(sci_not),'}');
        if offset ~= 0
            if sign(offset) > 0
                offset_sign = ' +';
            else
                offset_sign = ' -';
            end
            legend_key = strcat(legend_key, offset_sign, num2str(abs(offset)));
        end
        legend_keys(end+1) = legend_key;
      
    end
    
    set(gca, 'NextPlot', 'ReplaceChildren');
    set(gca, 'LineStyleOrder', {'-',':','--'});
    plot(include-1, all_props);
    set(gca, 'ytick', [-20:0.5:20]);
    legend(legend_keys);
  else
    disp('No extended properties found.');
    fig = [];
  end
end

function sci_not = sci_notation(value)
    if (value > 0)
        sci_not = ceil(log10(value)) - 1;
    else
        sci_not = 0;
    end
end

  