function add_vox_lines_to_plot(vox_lengths, num_voxels, colourize, offset, transparency)
% function add_vox_lines_to_plot(vox_size)
% Adds voxel markers to plot

  if length(vox_lengths) == 1
      vox_lengths = [vox_lengths, vox_lengths, vox_lengths];
  end

  if ~exist('colourize', 'var')
    colourize = true;
  end
  
  if ~exist('offset', 'var')
    offset = -half_vox_lengths.*num_voxels;
  end

  if ~exist('transparency', 'var')
    transparency = 1;
  end
  
  label_contrast = 0.8;
  base_colour = [.75,.75,.75];
  
  image_extent = vox_lengths .* num_voxels + offset;
  
  hold on;
  
  for dim_i = 1:3
      row_dim = mod(dim_i, 3) + 1;
      col_dim = mod(dim_i + 1, 3) + 1;
      for row_pos = offset(row_dim):vox_lengths(row_dim):image_extent(row_dim)
          for col_pos = offset(col_dim):vox_lengths(col_dim):image_extent(col_dim)
              
              if colourize && (row_pos == 0) && (col_pos == 0)
                  colour = [0.5,0.5,0.5];
                  colour(dim_i) = label_contrast;
              else
                  colour = base_colour;
              end

              line = [offset; image_extent];
              line(:, row_dim) = row_pos;
              line(:, col_dim) = col_pos;

              if transparency == 1
                  plot3(line(:, 1), line(:, 2), line(:, 3), 'Color', colour);
              else
                  patchline(line(:, 1), line(:, 2), line(:, 3), ...
                      'edgecolor', colour, 'edgealpha', transparency);
              end
              
          end
      end
  end

  hold off;   
  
end