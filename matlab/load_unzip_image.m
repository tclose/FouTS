function [image_matrix, labels] = load_unzip_image(filename)

  image = read_image(filename);
  
  image_matrix = image.data(:)';
  labels{1} = 'All Voxel Directions';
  
end