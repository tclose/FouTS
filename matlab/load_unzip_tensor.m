function [tensor_matrix, labels] = load_unzip_tensor(filename)
  
  [tensor_array, row_labels] = load_tensor(filename);

  [tensor_matrix, labels] = unzip_tensor(tensor_array, row_labels);
    
end