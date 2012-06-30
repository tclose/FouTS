function img = plot_neighbourhood(image_filename, mid_pos, figure_start, plot_limits)

	load('/data/home/tclose/Code/Tractography/matlab_scripts/grad_scheme');
	load('/data/home/tclose/Code/Tractography/matlab_scripts/p_scheme');

	if (~exist('image_filename'))
		image_filename = '~/Data/Set_of_fifty/chunked_separate_images/001-11-11-11.mif';
	end
	
	if (~exist('mid_pos'))
		mid_pos = [1 1 1];
	end
	
	if (~exist('figure_start'))
		figure_start = 0;
	else
		figure_start = figure_start - 1;
	end

	if (~exist('plot_limits'))
		dynamically_set_plot_limits = 1;
	else
		dynamically_set_plot_limits = 0;		
	end	
	

	img_struct = read_image(image_filename);

    
	
	if (~isfield(img_struct, 'data'))
		error(['Could not read image from file ' image_filename]);
	else

		img = img_struct.data;
		img = img(mid_pos(1):(mid_pos(1)+2), mid_pos(2):(mid_pos(2)+2), mid_pos(3):(mid_pos(3)+2), :);
		
		size_of_image = size(img);
		if (size_of_image(1:3) ~= [3 3 3])
			error('plot_neighbourhood currently only plots 3x3 neighbourhoods of voxels');
		end
		
		if (dynamically_set_plot_limits)
			plot_limits = max(max(max(max(abs(img)))));
    end

    if plot_limits == 0
      error('Neighbourhood contains no signal');
    end
    
		for (z = 1:3)
      
      next_fig = z+figure_start;
      if (ishandle(next_fig))
        
        close(next_fig);
      end
      
			fig = figure(next_fig);

			set(fig,'Units','normalized') 

			set(fig, 'Position', [0.05 + ((z-1) * 0.3), (0.55 * mod(figure_start+1,2) + 0.05) 0.3 0.375]);
			set(fig, 'DoubleBuffer', 'on');
			set(fig, 'Name', 'Amplitude plot');

			cameratoolbar('Show');
			cameratoolbar('SetMode','orbit');

			set(fig, 'Color', [1 1 1])

			for (y = 1:3)

				for (x = 1:3)

					
					ax = subplot(3,3, (y-1)*3 + x);
					
% 					title(['z=' num2str(z) ', y=' num2str(y) ', x=' num2str(x)]);
					
					set(ax, 'Xlim', [-plot_limits plot_limits]);
					set(ax, 'Ylim', [-plot_limits plot_limits]);
					set(ax, 'Zlim', [-plot_limits plot_limits]);
					
					voxel_signal = squeeze(img(x,y,z,:));
					sh = amp2SH(voxel_signal, grad_scheme);

					plot_SH(sh, p_scheme, ax);

					daspect([1,1,1]);
				end
			end
		end
	end
end


