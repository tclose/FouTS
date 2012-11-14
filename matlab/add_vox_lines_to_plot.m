function add_vox_lines_to_plot(vox_size, num_voxels, colourize, offset)
% function plot_vox_lines(vox_size)

  if ~exist('colourize')
    colourize = true;
  end

  label_contrast = 0.8;

  half_vox_size = vox_size/2;

  hold on;
  
  start_edge = -half_vox_size*num_voxels;
  end_edge   = half_vox_size*num_voxels;
  
  base_colour = [.75,.75,.75];
  
  for row_pos = start_edge:vox_size:end_edge
    for col_pos = start_edge:vox_size:end_edge
      
      line_mat = [row_pos, col_pos, start_edge; row_pos, col_pos, end_edge];
      line_mat(:,1) = offset + line_mat(:,1);
      line_mat(:,2) = offset + line_mat(:,2);     
      
      for dim_i = 0:2

        colour = base_colour;        
        
        if colourize && (row_pos == start_edge) && (col_pos == start_edge)
          colour = [0.5,0.5,0.5];
          colour(dim_i+1) = label_contrast;
        end
          
        plot3(line_mat(:,mod(0+dim_i,3)+1),line_mat(:,mod(1+dim_i,3)+1),line_mat(:,mod(2+dim_i,3)+1), 'Color', colour);
        
  
      end
        
    end
  end  

  hold off;   
  
end