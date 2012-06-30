function [tensor_matrix, labels] = unzip_tensor(tensor_array, row_labels)
  
  num_rows = size(tensor_array,1);
  
  if num_rows ~= size(tensor_array,2)
    error(['Number columns (' num2str(size(tensor_array,2)) ') does not match number of rows (' num2str(num_rows) ').']);
  end
  
  if num_rows ~= length(row_labels)
    error(['Length of column labels (' num2str(length(row_labels)) ') does not match number of rows (' num2str(num_rows) ').']);
  end
  
  labels = cell(num_rows^2,1);
  
  for row_i = 1:num_rows
    for col_i = 1:num_rows 
      
      row_label = strrep(row_labels(row_i), '_', ' ');
      col_label = strrep(row_labels(col_i), '_', ' ');
      
      row_label = row_label{1};
      col_label = col_label{1};
      
      row_label = row_label(1:(end-1));
      col_label = col_label(1:(end-1));
      
      labels{(row_i-1)*num_rows + col_i} = [row_label ' :: ' col_label];
    end
  end
  
  tensor_matrix = reshape(tensor_array, num_rows^2, size(tensor_array,3));
  

end