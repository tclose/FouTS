function [strand_matrix, labels] = unzip_strands(strands, prop_keys, prop_values, label_prefix)

  if ~exist('label_prefix', 'var')
    label_prefix = '';
  end

  num_strands = size(strands,1);

  if num_strands == 0
    error ('No strands loaded.');
  end

  degree = size(strands{1},1);

  for strand_i = 1:num_strands
    if size(strands{strand_i},1) ~= degree
      error(['Degree of strand_matrix strand ' num2str(strand_i) ' (' num2str(size(strands{strand_i},1)) ') does match that of previous strands  (' num2str(degree) ').']);
    end
  end

  num_props = length(prop_keys);    
  num_parameters = degree * 3 + num_props;

  strand_matrix = zeros(num_parameters, num_strands);
  labels = cell(num_parameters,1);

  for prop_i = 1:num_props
    strand_matrix(prop_i,:) = get_properties(prop_keys, prop_values, prop_keys{prop_i});
    labels{prop_i} = [label_prefix prop_keys{prop_i}];
  end

  for strand_i = 1:num_strands    
    
    start = num_props + 1;
%     finish = num_props + degree * 3;

    transposed_strand = strands{strand_i}';
    strand_matrix(start:end,strand_i) = transposed_strand(:);
    
  end

  for degree_i = 1:degree
    for dim_i = 1:3
      param_i = (degree_i-1) * 3 + dim_i + num_props;
      labels{param_i} = [label_prefix 'Degree ' num2str(degree_i-1) ' - Dim ' num2str(dim_i-1) ];
    end
  end
  
end