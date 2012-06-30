function fig = plot_extend_elem_properties(figure_name, all_prop_keys, all_prop_values, bundle_colours, include, elem_include, fig_index, num_figures, omit_properties)
%Not implemented yet.
       
  if ~exist('include','var')
    include = [];
  end
    
  if ~exist('elem_include','var')
    elem_include = [];
  end

  if ~exist('fig_index', 'var') 
    fig_index = 3;
  end
   
  if ~exist('num_figures', 'var') 
    num_figures = 3;
  end

  if ~exist('omit_properties', 'var')
    omit_properties = cell(0);
  end
  
  if ~isempty(all_prop_keys)

    if ~isempty(omit_properties)
      prop_keys = [];
      
      for key_i = 1:length(all_prop_keys)
        
        omit = 0;
        
        for omit_i = 1:length(omit_properties)
          if strfind(omit_properties{omit_i}, all_prop_keys{key_i})
            omit = 1;
            break;
          end
        end
        
        if ~omit
          prop_keys{length(prop_keys)+1} = all_prop_keys{key_i}; %#ok<AGROW>
        end
        
      end
    else
      prop_keys = all_prop_keys;
    end
      
    if isempty(include)
      include = [1:size(all_prop_values,1)];
    end

    
    fig = my_figure(figure_name,fig_index,num_figures,[]);

    num_keys = length(prop_keys);

    for key_i = 1:num_keys

      key = prop_keys(key_i);
      
      all_props = [];

      for set_i = include

        prop_values = all_prop_values{set_i};

        if ~isempty(elem_include)
          prop_values = prop_values(elem_include+1,:);
        end

        all_props = [all_props; get_properties(prop_keys, prop_values, key)'];

      end

      if strcmp(key,'tract_volume')
        all_props = all_props * (pi / 0.15^2);
        sum(all_props) ./ size(all_props,1)
      end
      
      num_strands = size(all_props,2);

      subplot(num_keys,1,key_i);
      set(gca, 'NextPlot', 'ReplaceChildren');

      if isempty(elem_include)
        colour_order = 1:num_strands;
      else
        colour_order = elem_include + 1;
      end

      set(gca, 'ColorOrder', bundle_colours(colour_order,:));
      plot(include, all_props);
      ylabel(strrep(key,'_',' '));

    end
        
  else
    disp('No extended, extended properties found');
    fig = [];
  end
end
