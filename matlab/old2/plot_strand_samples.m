function [samples, colours_of_strands] = plot_strand_samples(varargin)
%  
% PURPOSE: Plots strands from strand files, and optionally displays reference sphere and voxels
%  
% ARGUMENTS: 
%  
%           filename     The filename of the strands in either .tck or .frr formats.
%  
% OPTIONS (name, description, type, default):
%  
%          -transparency  
%                 transparency of the plotted samples
%                 float
%                 0.75
%  
%          -include        
%                 The indices of the strands to include in the plot
%                 matrix_1x:
%  
%          -sample_include 
%                 The indices of the samples to include in the plot
%                 matrix_1x:
%  
%          -colours_of_strands
%                 Colours of the plotted strands
%                 matrix_:x3
%  
%          -voxel_size     
%                 Size of reference voxel
%                 float
%                 0.15
%  
%          -num_voxels     
%                 Number of reference voxels
%                 int
%                 3
%  
%          -cube_size      
%                 Size of reference voxel
%                 float
%                 0
%  
%          -num_segments   
%                 Number of segments to plot for each strand.
%                 int
%                 100
%  
%          -tube_radius    
%                 Radius of tubes to plot for each strand.
%                 float
%                 0
%                 
                
  global colours_of_strands;

% The following code sets the input arguments and options as specified by the 'arguments' and 'options' cell arrays. 
  
  description = 'Plots strands from sample files, and optionally displays reference sphere and voxels';
  
  arguments = {'filename', 'The filename of the strand samples (''.smp'' format).'};
  
  options = {...
            'transparency  ', 0.75,          'float',  'transparency of the plotted samples';...
            'include        ', [],          'matrix_1x:', 'The indices of the strands to include in the plot';...            
            'sample_include ', [],          'matrix_1x:', 'The indices of the samples to include in the plot';...            
            'colours_of_strands', colours_of_strands,     'matrix_:x3', 'Colours of the plotted strands';... 
            'voxel_size     ', 0.15,        'float',  'Size of reference voxel';...
            'num_voxels     ', 3,           'int',    'Number of reference voxels';...
            'cube_size      ', 0,           'float',  'Size of reference voxel';...
            'num_segments   ', 100,         'int',    'Number of segments to plot for each strand.';...
            'tube_radius    ', 0,           'float',  'Radius of tubes to plot for each strand.';...
            ...%'sphere_radius  ', 0.389711,    'float',  'Size of reference sphere'...
            };         


  num_args = nargin;
  parse_arguments      
  if (help_display) 
    return;
  end

%   End arguments %


  
  dots_i = findstr(filename,'.');
	ext = filename(dots_i(end):end);
  
  if ~strcmp(ext, '.smp')
    error(['Input file must be a samples file (extension ''.smp''), found file with extension ' ext '.']);    
  end
  
  [samples, ext_prop_keys, ext_prop_values, initial_set, true_set, sphere_radius] = load_strand_samples(filename); 

  num_strands = size(initial_set,1);
  num_samples = size(samples,1) - 1;
  
  if any(find(sample_include > num_samples))
    error('Selected sample indices exceed number of loaded samples');
  end
  
  if isempty(sample_include)
    sample_include = 1:num_samples;
  end
  

  bundle_indices_index = find(strcmp(ext_prop_keys, 'bundle_index'));
  
  cdata = cell(num_samples, num_strands);
  
  
  if ~isempty(bundle_indices_index)

%     bundle_indices_strings = ext_prop_values{bundle_indices_index};
% 
%     bundle_indices = [];
%     
%     for sample_i = 1:size(bundle_indices_strings,1)
%       bundle_indices = [bundle_indices; regexp(bundle_indices_strings,'|','split')]; %#ok<AGROW>
%     end  
      
    bundle_indices = ext_prop_values;

    max_bundle_index = max(max(bundle_indices))+1;
    
    if size(colours_of_strands,1) < max_bundle_index %#ok<NODEF>
      colours_of_strands = rand([max_bundle_index, 3]); %#ok<NASGU>
      display_colour_key(colours_of_strands, 0:1:(max_bundle_index-1)); %#ok<NODEF>
    end
    
    for sample_i = 1:num_samples
    
      for strand_i = 1:num_strands

        colour_data = repmat(shiftdim(colours_of_strands(bundle_indices(sample_i, strand_i)+1,:),-1), [num_segments 1 1]); %#ok<NODEF>
        cdata{sample_i, strand_i} = [colour_data colour_data];

      end
      
    end
      
  else
    
    if size(colours_of_strands,1) < num_strands %#ok<NODEF>
      colours_of_strands = rand([num_strands, 3]); %#ok<NASGU>
      display_colour_key(colours_of_strands, 0:1:(num_strands-1)); %#ok<NODEF>
    end
    for sample_i = 1:num_samples
      for strand_i = 1:num_strands

        colour_data = repmat(shiftdim(colours_of_strands(strand_i,:),-1), [num_segments 1 1]); %#ok<NODEF>
        cdata{sample_i, strand_i} = [colour_data colour_data];

      end
    end
    
  end

  alpha_data = ones(num_segments,2) * transparency; 
  
	main_fig = figure();

	set(main_fig,'Units','normalized') 

	set(main_fig, 'Position', [0.3 0.25 0.4 0.5]);
	set(main_fig, 'DoubleBuffer', 'on');
	set(main_fig, 'Name', strrep(filename, '_', ' '));
	
	cameratoolbar('Show');
	cameratoolbar('SetMode','orbit');
  campos([0 0 1]);
  camup([0 1 0]);

	clf;

  whitebg(main_fig,'black');
	hold on;

  
  fprintf('\n');
	disp(['Plotting ' num2str(size(sample_include,2)) ' samples.']);
  
  for sample_i = sample_include

    sample = samples{sample_i,1};
    
%     if tube_radius > 0
%       
%       display_strands(sample, colours_of_strands, include, 'strand', '3D');      
%       
%     else
    
      for strand_i = 1:num_strands

        tck = fourier2tck(sample{strand_i,1});

        h = surface(...
          'XData', [tck(:,1) tck(:,1)],...
          'YData', [tck(:,2) tck(:,2)],...
          'ZData', [tck(:,3) tck(:,3)],...
          'CData', cdata{sample_i, strand_i},...
          'AlphaData', alpha_data,...
          'FaceColor', 'none',...
          'EdgeColor', 'flat',...
          'Marker', 'none', ...
          'LineStyle', '-', ...        
          'LineWidth', 1);
        view(3)     
        daspect ([ 1 1 1 ]);
      end  

      if ~mod(sample_i,25)
        fprintf('.');
      end

%     end
      
      
      
      
  end
  
	hold off;

	%set(gca, 'color', [0 0 0]);
	set(gcf, 'color', [0 0 0]);
	
	label_handle = get(gca,'xlabel');
	set(label_handle, 'string', 'X-axis', 'color', [1 1 1]);

  label_handle = get(gca,'ylabel');
	set(label_handle, 'string', 'Y-axis', 'color', [1 1 1]);
   
	label_handle = get(gca,'zlabel');
	set(label_handle, 'string', 'Z-axis', 'color', [1 1 1]);
	
	h = get (gca, 'children');
	daspect ([ 1 1 1 ]);

 	light           
 	lighting gouraud;
	

%	figure(main_fig);


  
%   if sphere_radius ~= 0
%     display_reference_sphere(sphere_radius, 0.1);
% 	end	
  
  if voxel_size ~= 0
    plot_vox_lines(voxel_size,num_voxels);
  end
  

  samples = samples(sample_include);
  
end