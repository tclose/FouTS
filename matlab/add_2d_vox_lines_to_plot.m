function add_2d_vox_lines_to_plot(vox_size, num_voxels, colour, offset)
% function plot_vox_lines(vox_size)

  if ~exist('vox_size','var')
    vox_size = 0.15;
  end

  if ~exist('num_voxels','var')
    num_voxels = 3;
  end
  
  if ~exist('colour','var')
    colour = 'white';
  end

  if ~exist('offset','var')
    offset = 0;
  end
  
  half_vox_size = vox_size/2;

  hold on;
  
  start_edge = -half_vox_size*num_voxels;
  end_edge   = half_vox_size*num_voxels;
  
  for row_pos = start_edge:vox_size:end_edge
    for col_pos = start_edge:vox_size:end_edge
      
      line_mat = [row_pos, col_pos, start_edge; row_pos, col_pos, end_edge] + offset;
      
      for dim_i = 0:1
          
        plot(line_mat(:,mod(1+dim_i,3)+1),line_mat(:,mod(2+dim_i,3)+1), 'Color', colour);
        
      end
        
    end
  end  

  hold off;   
  
end