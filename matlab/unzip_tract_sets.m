function [tract_sets_matrix, labels] = unzip_tract_sets(sets, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values)

  num_sets = size(sets,1);

  if num_sets == 0
    error ('No sets loaded.');
  end

  num_tracts = size(sets{1},1);
  first_set = sets(1);
  degree = size(first_set{1}{1},1);

  for set_i = 1:num_sets
    
    set = sets{set_i};
    new_num_tracts = size(set,1);
    
    if new_num_tracts ~= num_tracts
      error(['Number of tracts in set ' num2str(set_i) ' (' num2str(new_num_tracts) ') does match that of previous sets  (' num2str(num_tracts) ').']);
    end
     
    for tract_i = 1:num_tracts
      
      new_degree = size(set{tract_i,1},1);
      
      if new_degree ~= degree
        error(['Degree of tract ' num2str(tract_i) ' in set ' num2str(set_i) ' (' num2str(new_degree) ') does match that of previous tracts  (' num2str(degree) ').']);
      end
    end
     
  end

  
  num_set_props = length(set_prop_keys);
  num_elem_props = length(elem_prop_keys);
  
  num_parameters = num_tracts * degree * 3 + num_set_props + num_elem_props * num_tracts;


  tract_sets_matrix = zeros(num_parameters, num_sets);
  labels = cell(num_parameters,1);

  for prop_i = 1:num_set_props
    tract_sets_matrix(prop_i,:) = get_properties(set_prop_keys, set_prop_values, set_prop_keys{prop_i});
    labels{prop_i} = set_prop_keys{prop_i};
  end

  
  for tract_i = 1:num_tracts
    
    same_tract_across_sets = cell(num_sets,3);
    same_tract_props_across_sets = cell(num_sets, num_elem_props);
    
    for set_i = 1:num_sets
      
      set = sets{set_i};
      
      same_tract_across_sets(set_i,:) = set(tract_i,:);
      same_tract_props_across_sets(set_i,:) = elem_prop_values{set_i}(tract_i,:);
    end
    
    [tract_matrix, tract_labels] = unzip_tracts(same_tract_across_sets, elem_prop_keys, same_tract_props_across_sets, ['Tract ' num2str(tract_i-1) ' - ']);
    
    start = (tract_i-1) * (degree * 9 + num_elem_props) + num_set_props + 1;
    finish = tract_i * (degree * 9 + num_elem_props) + num_set_props;
    
    tract_sets_matrix(start:finish,:) = tract_matrix;
    labels(start:finish) = tract_labels;
    
  end
  
end