function main_fig = plot_tcks(varargin)
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
%          -colours_of_bundles
%                 Colours of the plotted tcks
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

  global colours_of_bundles;

  description = 'Plots tcks from tcks files, and optionally displays reference sphere and voxels';
  
  arguments = {'tcks_filename', 'The filename of the tcks in either .tck or .frr formats.'};
  
  options = {...
            'tcks_radius  ', 0.02,        'float',  'Radii of the plotted tcks';...
            'num_points     ', 100,         'int',    'Number of points to plot along each tcks.';...            
            'include        ', [],          'matrix_:x:', 'The indices of the tcks to include in the plot';...            
            'bundle_include ', [],          'matrix_:x:', 'The indices of the bundles to include in the plot';...            
            'colours_of_bundles',        colours_of_bundles,     'matrix_:x3', 'Colours of the plotted tcks';...
            'voxel_size     ', 0.15,        'float',  'Size of reference voxel';...
            'num_voxels     ', 3,           'int',    'Number of reference voxels';...
            'cube_size      ', 0,           'float',  'Size of reference voxel';...
            'style          ', 'tubes',   'string', 'Plot to style (either ''tcks'' or ''lines'')..';...
            'line_style     ', '-',   'string', 'Line plot to style (internal Matlab linespec format). Only valid with ''-style lines''..';...
            'sphere_radius  ', 0,    'float',  'Size of reference sphere';...
            'no_axes',       0,    'bool', 'Removes axes from plot';...
            'tube_corners   ', 6,          'int',    'Change plot to ''line style''.(NB: Only relevant with ''-num_strands'' option = 0).';...
            'clean',         0,      'bool', 'Removes everything else from plot';...
            'happy_colours',         0,     'bool', 'Uses ''happy colours''';...
            'inv_happy_colours',         0,     'bool', 'Uses the inverse of ''happy colours''';...
            'no_voxline_highlight', 0, 'bool', 'Doesn''t highlight corner axes of voxel lines.';...
            'campos',      [],      'matrix_1x3', 'Position of the camera'};         


  parse_arguments      
  if (help_display) 
    return;
  end

  if clean
    voxel_size = 0;
    sphere_radius = 0;
    no_axes = 1;
  end
  
%   End arguments %
  
  if happy_colours
    
    if inv_happy_colours
      error('Can''t use ''-happy_colours'' and ''-inv_happy_colours'' simultaneously');
    end
    
    load('/home/tclose/Data/Tractography/misc/comb_happy_colours.mat');
    
  elseif inv_happy_colours
    
    load('/home/tclose/Data/Tractography/misc/inv_comb_happy_colours.mat');
    
  end

  [tcks, props, prop_keys, prop_values] = load_tcks(tcks_filename);
  
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
    
  add_colour_key(bundle_indices,colours_of_bundles);
  
  %Set up the figure.
  main_fig = my_figure(strrep(tcks_filename,'_',''));

  cameratoolbar('Show');
  cameratoolbar('SetMode','orbit');
  
  if strfind('tubes', style) == 1
    
    add_tcks_to_plot(tcks, radii, colours_of_bundles, bundle_indices, tube_corners); 
    
  elseif strfind('lines', style) == 1
    
    add_lines_to_plot(tcks, colours_of_bundles, bundle_indices, line_style);    
    
  else
    
    error(['Urecogised style option ''' style '''.']);
    
  end
    
  
  if sphere_radius ~= 0
    
    add_sphere_to_plot(sphere_radius);
    
	end	
  
  if voxel_size ~= 0
    
    add_vox_lines_to_plot(voxel_size,num_voxels,~no_voxline_highlight);
    
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