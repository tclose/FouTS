function [tract_matrix, labels] = unzip_tracts(tracts, prop_keys, prop_values, label_prefix)

  if ~exist('label_prefix', 'var')
    label_prefix = '';
  end

  num_tracts = size(tracts,1);

  if num_tracts == 0
    error ('No tracts loaded.');
  end

  degree = size(tracts{1,1},1);

  for tract_i = 1:num_tracts
    for ax_i = 1:3
      if size(tracts{tract_i,ax_i},1) ~= degree
        error(['Degree of tract_matrix tract ' num2str(tract_i) ' (' num2str(size(tracts{tract_i,ax_i},1)) ') does match that of previous tracts  (' num2str(degree) ').']);
      end
    end
  end

  num_props = length(prop_keys);    
  num_parameters = degree * 3 * 3 + num_props;

  tract_matrix = zeros(num_parameters, num_tracts);
  labels = cell(num_parameters,1);

  for prop_i = 1:num_props
    tract_matrix(prop_i,:) = get_properties(prop_keys, prop_values, prop_keys{prop_i});
    labels{prop_i} = [label_prefix prop_keys{prop_i}];
  end

  for tract_i = 1:num_tracts   
    for ax_i = 1:3        
      start = num_props + (ax_i-1) * degree * 3 + 1;
      finish = num_props + (ax_i) * degree * 3;
      transposed_tract = tracts{tract_i,ax_i}';
      tract_matrix(start:finish,tract_i) = transposed_tract(:);
    end
  end

  for ax_i = 0:1:2
    for degree_i = 0:1:(degree-1)
      for dim_i = 0:1:2
        param_i = ax_i*degree*3 + degree_i * 3 + dim_i + (num_props+1);
        labels{param_i} = [label_prefix 'Axis ' num2str(ax_i) ' - Degree ' num2str(degree_i) ' - Dim ' num2str(dim_i) ];
      end
    end
  end

  
end