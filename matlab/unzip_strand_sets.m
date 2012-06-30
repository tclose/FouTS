function [strand_sets_matrix, labels] = unzip_strand_sets(sets, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values)

  num_sets = size(sets,1);

  if num_sets == 0
    error ('No sets loaded.');
  end

  num_strands = size(sets{1},1);
  first_set = sets(1);
  degree = size(first_set{1}{1},1);

  for set_i = 1:num_sets
    
    set = sets{set_i};
    
    if size(sets{set_i},1) ~= num_strands
      error(['Number of strands in set ' num2str(set_i) ' (' num2str(size(sets{set_i},1)) ') does match that of previous sets  (' num2str(num_strands) ').']);
    end
     
    for strand_i = 1:num_strands
      if size(sets{set_i}{strand_i},1) ~= degree
        error(['Degree of strand ' num2str(strand_i) ' in set ' num2str(set_i) ' (' num2str(size(sets{set_i}{strand_i},1)) ') does match that of previous strands  (' num2str(degree) ').']);
      end
    end
     
  end

  
  num_set_props = length(set_prop_keys);
  num_elem_props = length(elem_prop_keys);
  
  num_parameters = num_strands * degree * 3 + num_set_props + num_elem_props * num_strands;


  strand_sets_matrix = zeros(num_parameters, num_sets);
  labels = cell(num_parameters,1);

  for prop_i = 1:num_set_props
    strand_sets_matrix(prop_i,:) = get_properties(set_prop_keys, set_prop_values, set_prop_keys{prop_i});
    labels{prop_i} = set_prop_keys{prop_i};
  end

  
  for strand_i = 1:num_strands
    
    same_strand_across_sets = cell(num_sets);
    same_strand_props_across_sets = cell(num_sets, num_elem_props);
    
    for set_i = 1:num_sets
      
      set = sets{set_i};
      
      same_strand_across_sets(set_i) = set(strand_i);
      
      if num_elem_props
        same_strand_props_across_sets(set_i) = elem_prop_values{set_i}(strand_i,:);
      end
    end
    
    [strand_matrix, strand_labels] = unzip_strands(same_strand_across_sets, elem_prop_keys, same_strand_props_across_sets, ['Strand ' num2str(strand_i-1) ' - ']);
    
    start = (strand_i-1) * (degree * 3 + num_elem_props) + num_set_props + 1;
    finish = strand_i * (degree * 3 + num_elem_props) + num_set_props;
    
    strand_sets_matrix(start:finish,:) = strand_matrix;
    labels(start:finish) = strand_labels;
    
  end
  
end