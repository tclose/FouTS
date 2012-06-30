function main_fig = plot_strands(varargin)
%  
% PURPOSE: Plots strands from strand files, and optionally displays reference sphere and voxels
%  
% ARGUMENTS: 
%  
%           strands_filename     The filename of the strands in either .tck or .frr formats.
%  
% OPTIONS (name, description, type, default):
%  
%          -strand_radius  
%                 Radii of the plotted strands
%                 float
%                 0.02
%  
%          -num_points     
%                 Number of points to plot along each strand.
%                 int
%                 100
%  
%          -include        
%                 The indices of the strands to include in the plot
%                 matrix_:x:
%  
%          -bundle_include 
%                 The indices of the bundles to include in the plot
%                 matrix_:x:
%  
%          -colours_of_bundles
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
%          -style          
%                 Plot to style (either 'strands' or 'lines')..
%                 string
%                 'strands'
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
%                 matrix_1:3

  main_fig = -1;

  global colours_of_bundles;

  description = 'Plots strands from strand files, and optionally displays reference sphere and voxels';
  
  arguments = {'strands_filename', 'The filename of the strands in either .tck or .frr formats.'};
  
  options = {...
            'strand_radius  ', 0.02,        'float',  'Radii of the plotted strands';...
            'num_points     ', 100,         'int',    'Number of points to plot along each strand.';...            
            'include        ', [],          'matrix_:x:', 'The indices of the strands to include in the plot';...            
            'bundle_include ', [],          'matrix_:x:', 'The indices of the bundles to include in the plot';...            
            'colours_of_bundles', colours_of_bundles,     'matrix_:x3', 'Colours of the plotted strands';...
            'voxel_size     ', 0.15,        'float',  'Size of reference voxel';...
            'num_voxels     ', 3,           'int',    'Number of reference voxels';...
            'cube_size      ', 0,           'float',  'Size of reference voxel';...
            'style          ', 'tubes',   'string', 'Plot to style (either ''strands'' or ''lines'')..';...
            'sphere_radius  ', 0,    'float',  'Size of reference sphere';...
            'sphere_transparency', 0.25,    'float',  'Transparency of reference sphere';...
            'tube_corners   ', 6,          'int',    'Change plot to ''line style''.(NB: Only relevant with ''-num_strands'' option = 0).';...
            'no_axes',       0,    'bool', 'Removes axes from plot';...
            'happy_colours',         0,     'bool', 'Uses ''happy colours''';...
            'inv_happy_colours',         0,     'bool', 'Uses the inverse of ''happy colours''';...
            'clean',         0,      'bool', 'Removes everything else from plot';...
            'campos',        [],      'matrix_1:3', 'Position of the camera';...
            'no_voxline_highlight', 0, 'bool', 'Doesn''t highlight corner axes of voxel lines.';...
            'invisible', 0,     'bool', 'Makes the figure invisible (for automatically saving afterwards).'};         


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


  [strands, props, prop_keys, prop_values] = load_strands(strands_filename);
  
   
  if ~isempty(bundle_include)
  
    num_strands = size(strands,1);
    
    if ~isempty(include)
      error('Cannot use ''include'' and ''bundle_include simultaneously');
    end    
    
    bundle_indices  = get_properties(prop_keys, prop_values, 'bundle_index', 0:1:(num_strands-1), num_strands);
    
    for strand_i = 1:num_strands
      if any(find(bundle_include == bundle_indices(strand_i)))
        include = [include; strand_i]; %#ok<AGROW>
      end
    end
    
    if isempty(include)
      disp('Warning! No bundles matched the given indices, returning all strands.');
    end
    
    %To match 0 indexing
    include = include - 1;
    
  end
  
  
  if ~isempty(include)
    
    strands = strands(include+1);
    prop_values = prop_values(include+1,:);
  
  end
  
  num_strands = size(strands,1);

  bundle_indices  = get_properties(prop_keys, prop_values, 'bundle_index', [0:1:(num_strands-1)]', num_strands);
  radii           = get_properties(prop_keys, prop_values, 'track_radius', strand_radius, num_strands);
    
  add_colour_key(bundle_indices,colours_of_bundles);
  
  main_fig = my_figure(strrep(strands_filename,'_',''), 1, 1, [1 1 1], 1, [],[],~invisible);

  cameratoolbar('Show');
  cameratoolbar('SetMode','orbit');
  
  
  if strcmp(file_extension(strands_filename), 'str')
    tcks = strands2tcks(strands);
  elseif strcmp(file_extension(strands_filename), 'tck')
    tcks = strands;
  else
    error (['Unrecognised extension, ''' extension(strands_filename) ''', of input file.']);
  end
  
  if strfind('tubes', style) == 1
    
    add_tcks_to_plot(tcks, radii, colours_of_bundles, bundle_indices, tube_corners); 
    
  elseif strfind('lines', style) == 1
    
    add_lines_to_plot(tcks, colours_of_bundles, bundle_indices);    
    
  else
    
    error(['Urecogised style option ''' style '''.']);
    
  end
    
  
  if sphere_radius ~= 0  
    
    add_sphere_to_plot(sphere_radius, sphere_transparency);
    
	end	
  
  if voxel_size ~= 0
    
    add_vox_lines_to_plot(voxel_size,num_voxels,~no_voxline_highlight);
    
  end
  
  
  set(gcf, 'color', [0 0 0]);

  if no_axes
    
    set(gca, 'visible', 'off') ;
    
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
  disp(['Plotted ' num2str(num_strands) ' strands.']);   
  
end