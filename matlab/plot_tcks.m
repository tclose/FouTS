function [main_fig, colour_indices] = plot_tcks(varargin)
%  
% PURPOSE: Plots tcks from tcks files, and optionally displays reference sphere and voxels
%  
% ARGUMENTS: 
%  
%           tcks_filename     The filename of the tcks in either .tck or .frr formats.
%  
% OPTIONS (name, description, type, default):
%  
%          -tcks_radius  
%                 Radii of the plotted tcks
%                 float
%                 0.02
%  
%          -num_points     
%                 Number of points to plot along each tcks.
%                 int
%                 100
%  
%          -include        
%                 The indices of the tcks to include in the plot
%                 matrix_:x:
%  
%          -bundle_include 
%                 The indices of the bundles to include in the plot
%                 matrix_:x:
%  
%          -colours
%                 Colours of the plotted tcks
%                 matrix_:x3
%  
%          -voxel_lengths     
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
%          -style          
%                 Plot to style (either 'tcks' or 'lines')..
%                 string
%                 'tcks'
%  
%          -sphere_radius  
%                 Size of reference sphere
%                 float
%                 0
%  
%          -no_axes
%                 Removes axes from plot
%                 bool
%                 0
%  
%          -clean
%                 Removes everything else from plot
%                 bool
%                 0
%  
%          -campos
%                 Position of the camera
%                 matrix_1x3
% 

  main_fig = -1;
  description = 'Plots tcks from tcks files, and optionally displays reference sphere and voxels';
  
  arguments = {'tcks_filename', 'The filename of the tcks in either .tck or .frr formats.'};
  
  options = {...
            'tcks_radius  ', 0.02,        'float',  'Radii of the plotted tcks';...
            'num_points     ', 100,         'int',    'Number of points to plot along each tcks.';...            
            'include        ', [],          'matrix_:x:', 'The indices of the tcks to include in the plot';...            
            'bundle_include ', [],          'matrix_:x:', 'The indices of the bundles to include in the plot';...            
            'colours',        [],     'matrix_:x3', 'Colours of the plotted tcks';...
            'colour_indices',           [],         'matrix_:x1',   'The mapping from bundle indices to the colours used';...
            'compact_colours',  0,          'bool',   'Compact the range of colours created to only use the number needed';...
            'num_voxels     ',          [],         'matrix_1x3',   'Number of reference voxels';...
            'voxel_lengths',            [],         'matrix_1x3',   'Size of reference voxel';...
            'voxel_offsets  ',          [],         'matrix_1x3',   'Offset of reference voxels';...
            'voxel_transparency',       0.25,       'float',        'The alpha value of the voxel lines';...
            'style          ', 'tubes',   'string', 'Plot to style (either ''tcks'' or ''lines'')..';...
            'line_style     ', '-',   'string', 'Line plot to style (internal Matlab linespec format). Only valid with ''-style lines''..';...
            'sphere_radius  ', 0,    'float',  'Size of reference sphere';...
            'no_axes',       0,    'bool', 'Removes axes from plot';...
            'tube_corners   ', 6,          'int',    'Change plot to ''line style''.(NB: Only relevant with ''-num_strands'' option = 0).';...
            'clean',         0,      'bool', 'Removes everything else from plot';...
            'happy_colours',         0,     'bool', 'Uses ''happy colours''';...
            'inv_happy_colours',         0,     'bool', 'Uses the inverse of ''happy colours''';...
            'no_voxline_highlight', 0, 'bool', 'Doesn''t highlight corner axes of voxel lines.';...
            'obs_image',                [],         'string',       'Overlays a slice along the x observed image.';...
            'transparency'   ,          1,          'float',        'Fix the transparency of the tracts to a set value.  If 0 the intensities of the of the tracts are used instead.';...
            'campos',      [],      'matrix_1x3', 'Position of the camera'};         


  parse_arguments      
  if (help_display) 
    return;
  end

  if clean
    voxel_lengths = 0;
    sphere_radius = 0;
    no_axes = 1;
  end
  
%   End arguments %
  colours = get_happy_colours(colours, happy_colours, inv_happy_colours); %#ok<NODEF>
  
  [tcks, props, prop_keys, prop_values] = load_tcks(tcks_filename);
  
  if isempty(num_voxels)    %#ok<NODEF>
    [num_voxels, ~, voxel_lengths, voxel_offsets, obs_image] = get_observed_properties(props, obs_image); %#ok<NODEF,ASGLU>
    if ~isempty(num_voxels)
        num_voxels = num_voxels(1:3);
        voxel_lengths = voxel_lengths(1:3);
    end
  elseif isempty(voxel_offsets) %#ok<NODEF>
      voxel_offsets = - num_voxels .* voxel_lengths ./2;
  end
  
  if ~isempty(bundle_include)
  
    num_tcks = size(tcks,1);
    
    if ~isempty(include)
      error('Cannot use ''include'' and ''bundle_include simultaneously');
    end    
    
    bundle_indices  = get_properties(prop_keys, prop_values, 'bundle_index', 0:1:(num_tcks-1), num_tcks);    
    
    for tcks_i = 1:num_tcks
      if any(find(bundle_include == bundle_indices(tcks_i)))
        include = [include; tcks_i]; %#ok<AGROW>
      end
    end
    
    if isempty(include)
      disp('Warning! No bundles matched the given indices, returning all tcks.');
    end
  
    %To match 0 indexing
    include = include - 1;
    
  end
  
  
  if ~isempty(include)
    
    tcks = tcks(include+1);
    prop_values = prop_values(include+1,:);
  
  end
  
  num_tcks = size(tcks,1);

  bundle_indices  = get_properties(prop_keys, prop_values, 'bundle_index', [0:1:(num_tcks-1)]', num_tcks);
  radii           = get_properties(prop_keys, prop_values, 'track_radius', tcks_radius, num_tcks);
    
  [colours, colour_indices] = set_bundle_colours(bundle_indices, colours, colour_indices, compact_colours); %#ok<NODEF>
  colour_indices = add_colour_key(bundle_indices, colours, colour_indices); %#ok<NODEF>
  
  %Set up the figure.
  main_fig = my_figure(strrep(tcks_filename,'_',''));

  cameratoolbar('Show');
  cameratoolbar('SetMode','orbit');
  
  if strfind('tubes', style) == 1
    
    add_tcks_to_plot(tcks, radii, colours, bundle_indices, tube_corners, transparency, colour_indices); 
    
  elseif strfind('lines', style) == 1
    
    add_lines_to_plot(tcks, colours, bundle_indices, line_style, colour_indices);    
    
  else
    
    error(['Urecogised style option ''' style '''.']);
    
  end
    
  
  if sphere_radius ~= 0
    
    add_sphere_to_plot(sphere_radius);
    
	end	
  
  if ~isempty(num_voxels)

    add_vox_lines_to_plot(voxel_lengths,num_voxels,~no_voxline_highlight,...
          voxel_offsets, voxel_transparency);

  end
  
  
  set(gcf, 'color', [0 0 0]);


  if no_axes
    
    set(gca, 'visible', 'off');
    
  else

    label_handle = get(gca,'xlabel');
    set(label_handle, 'string', 'X-axis', 'color', [1 1 1]);

    label_handle = get(gca,'ylabel');
    set(label_handle, 'string', 'Y-axis', 'color', [1 1 1]);

    label_handle = get(gca,'zlabel');
    set(label_handle, 'string', 'Z-axis', 'color', [1 1 1]);
  
  end
  
  if ~isempty(campos)
    
    set(gca,'CameraPosition',campos);
    
  end
  
  daspect ([ 1 1 1 ]);

  light           
  lighting gouraud;


  fprintf('\n');
  disp(['Plotted ' num2str(num_tcks) ' tcks.']);   
  
end