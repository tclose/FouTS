function plot_vox_lines(vox_size, num_voxels)
% function plot_vox_lines(vox_size)

  label_contrast = 0.8;

  half_vox_size = vox_size/2;

  hold on;
  
  start_edge = -half_vox_size*num_voxels;
  end_edge   = half_vox_size*num_voxels;
  
  base_colour = [.75,.75,.75];
  
  for row_pos = start_edge:vox_size:end_edge
    for col_pos = start_edge:vox_size:end_edge
      
      line_mat = [row_pos, col_pos, start_edge; row_pos, col_pos, end_edge];
      
      for dim_i = 0:2

        colour = base_colour;        
        if (row_pos == start_edge) && (col_pos == start_edge)
          colour = [0.5,0.5,0.5];
          colour(dim_i+1) = label_contrast;
        end
          
        plot3(line_mat(:,mod(0+dim_i,3)+1),line_mat(:,mod(1+dim_i,3)+1),line_mat(:,mod(2+dim_i,3)+1), 'Color', colour);
        
  
      end
        
    end
  end  

  hold off;   
  
end