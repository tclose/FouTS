function [strand_sets, props, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values]  = load_strand_sets(filename, include, load_ext_props)

  if ~exist('include')
    include = [];
  end

  if ~exist('load_ext_props')
    load_ext_props = false;
  end

  [elems, props] = read_fibres(filename);
  
  sets = split_at_file_seperator(elems, [-inf,nan,inf]);
  
  num_sets = length(sets);
  
   
  if isempty(include)
    include = 1:num_sets;
    num_include = num_sets;
  else
    num_include = size(include(:),1);
  end
  
  strand_sets = cell(num_include,1);
  
  for set_i = include
    strand_sets{set_i} = parse_strands(sets{set_i});
  end
  
  if load_ext_props && exist([filename 'x'])
    [set_prop_keys, set_prop_values] = read_element_properties([filename 'x']);
    set_prop_values = set_prop_values(include',:);

  else
    
    set_prop_keys = cell(0);
    set_prop_values = cell(0);

    
  end 
  
  
  if load_ext_props && exist([filename 'xx'])
    [elem_prop_keys, elem_prop_values] = read_set_element_properties([filename 'xx']);    
    elem_prop_values = elem_prop_values(include',:);
    
  else
    
    elem_prop_keys = cell(0);
    elem_prop_values = cell(0);
    
  end  
  
end