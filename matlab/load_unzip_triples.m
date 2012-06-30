function [triples_matrix, labels] = load_unzip_triples(filename)

  triples_matrix = load_triples(filename)';
  
  labels = cell(3,1);
  
  labels{1} = 'Dim 0';
  labels{2} = 'Dim 1';
  labels{3} = 'Dim 2';
  
end