function main_fig = plot_tracts(varargin)
%  
% PURPOSE: Plots strands from strand files, and optionally displays reference sphere and voxels
%  
% ARGUMENTS: 
%  
%           tracts_filename     The filename of the tracts ('.tct' format).
%  
% OPTIONS (name, description, type, default):
%  
%          -include        
%                 The indices of the strands to include in the plot
%                 matrix_:x:
%  
%          -colours_of_bundles
%                 Colours of the plotted strands
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
%          -num_length_sections   
%                 Number of segments to plot for each tract.
%                 int
%                 100
%  
%          -strand_radius  
%                 Radii of the plotted strands. (NB: can only be used with '-style' option 'strands').
%                 float
%                 0.02
%  
%          -num_width_sections    
%                 Number of strands to plot along each axis. If zero the default 3D surface option will be used instead.
%                 int
%                 2
%  
%          -style          
%                 Change plot to 'line style'.(NB: Only relevant with '-num_width_sections' option > 0).
%                 string
%                 'tracts'
%  
%          -tube_corners   
%                 Change plot to 'line style'.(NB: Only relevant with '-num_width_sections' option = 0).
%                 int
%                 12
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
%          -colour_key
%                 Displays colour key
%                 bool
%                 0
%  
%          -oblong
%                 Prints an oblong tractlet rather than cylindrical
%                 bool
%                 0
%  
%          -happy_colours
%                 Uses 'happy colours'
%                 bool
%                 0
%  
%          -inv_happy_colours
%                 Uses the inverse of 'happy colours'
%                 bool
%                 0
%  
%          -ignore_bundles
%                 Ignores bundle indices when plotting colours.
%                 bool
%                 0
%  
%          -no_voxline_highlight
%                 Doesn't highlight corner axes of voxel lines.
%                 bool
%                 0
%  
%          -transparency
%                 Fix the transparency of the tracts to a set value.  If 0 the intensities of the of the tracts are used instead.
%                 float
%                 1
%  
%          -highlight_axes
%                 Highlights the axes when printing in 'strand' or 'line' style
%                 bool
%                 0

  main_fig = -1;

  global colours_of_bundles;

  description = 'Plots strands from strand files, and optionally displays reference sphere and voxels';
  
  arguments = {'tracts_filename', 'The filename of the tracts (''.tct'' format).'};
  
  options = {...
            'include        ',          [],         'matrix_:x:',   'The indices of the strands to include in the plot';...                       
            'colours_of_bundles',       colours_of_bundles,         'matrix_:x3', 'Colours of the plotted strands';...
            'num_voxels     ',          [],         'matrix_1x3',   'Number of reference voxels';...
            'voxel_lengths',            [],         'matrix_1x3',   'Size of reference voxel';...
            'voxel_offsets  ',          [],         'matrix_1x3',   'Offset of reference voxels';...
            'voxel_transparency',       0.25,       'float',        'The alpha value of the voxel lines';...
            'cube_size      ',          0,          'float',        'Size of reference voxel';...
            'num_length_sections   ',   100,        'int',          'Number of segments to plot for each tract.';...            
            'strand_radius  ',          0.02,       'float',        'Radii of the plotted strands. (NB: can only be used with ''-style'' option ''strands'').';...
            'num_width_sections    ',   2,          'int',          'Number of strands to plot along each axis. If zero the default 3D surface option will be used instead.';...            
            'strands_per_acs',          -1,         'float',        'Instead of a fixed number of width sections the number of plotted strands is determined by the ACS of the tract';...
            'style          ',          'tractlets','string',       'Change plot to ''line style''.(NB: Only relevant with ''-num_width_sections'' option > 0).';...
            'tube_corners   ',          12,         'int',          'Change plot to ''line style''.(NB: Only relevant with ''-num_width_sections'' option = 0).';...
            'sphere_radius  ',          0,          'float',        'Size of reference sphere';...
            'no_axes',                  0,          'bool',         'Removes axes from plot';...
            'clean',                    0,          'bool',         'Removes everything else from plot';...
            'colour_key',               0,          'bool',         'Displays colour key';...
            'oblong',                   0,          'bool',         'Prints an oblong tractlet rather than cylindrical';...
            'happy_colours',            0,          'bool',         'Uses ''happy colours''';...
            'inv_happy_colours',        0,          'bool',         'Uses the inverse of ''happy colours''';...
            'ignore_bundles',           0,          'bool',         'Ignores bundle indices when plotting colours.';...
            'no_voxline_highlight',     0,          'bool',         'Doesn''t highlight corner axes of voxel lines.';...
            'transparency'   ,          1,          'float',        'Fix the transparency of the tracts to a set value.  If 0 the intensities of the of the tracts are used instead.';...
            'highlight_axes',           0,          'bool',         'Highlights the axes when printing in ''strand'' or ''line'' style';...
            'invisible',                0,          'bool',         'Makes the figure invisible (for automatically saving afterwards).';...
            'hold_on'          ,        0,          'bool',         'plot on previous figure';...
            'obs_image',                [],         'string',       'Overlays a slice along the x observed image.';...
            'slice_x',                  [],         'matrix_:x:',   'Overlays slices of observed image along the given indices.';...
            'slice_y',                  [],         'matrix_:x:',   'Overlays slices of observed image along the given indices.';...
            'slice_z',                  [],         'matrix_:x:',   'Overlays slices of observed image along the given indices.'};         


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

  if happy_colours
    
    if inv_happy_colours
      error('Can''t use ''-happy_colours'' and ''-inv_happy_colours'' simultaneously');
    end
    
    load('/home/tclose/Data/Tractography/misc/comb_happy_colours.mat');
    
  elseif inv_happy_colours
    
    load('/home/tclose/Data/Tractography/misc/inv_comb_happy_colours.mat');
    
  end

  if ~strcmp(file_extension(tracts_filename), 'tct')
    error (['Extension, ''' file_extension(tracts_filename) ''' is not a valid tract file (''.tct'').']);
  end

  if (highlight_axes && ~mod(num_width_sections,2))
    error(['If ''-highlight_axes'' option is selected the number of width sections ''-num_width_sections'' needs to be odd (' num2str(num_width_sections) ' provided)'])
  end
  
  [tracts, props, prop_keys, prop_values] = load_tracts(tracts_filename);
   
  num_loaded_tracts = size(tracts,1);
  
  if isempty(num_voxels)    %#ok<NODEF>
    [num_voxels, ~, voxel_lengths, voxel_offsets, obs_image] = get_observed_properties(props, obs_image);
    if ~isempty(num_voxels)
        num_voxels = num_voxels(1:3);
        voxel_lengths = voxel_lengths(1:3);
    end
  elseif isempty(voxel_offsets) %#ok<NODEF>
      voxel_offsets = - num_voxels .* voxel_lengths ./2;
  end
  
  acs = get_properties(prop_keys, prop_values, 'alpha', (0:(num_loaded_tracts-1))', num_loaded_tracts).^2;  
  bundle_indices = get_properties(prop_keys, prop_values, 'bundle_index', (0:(num_loaded_tracts-1))', num_loaded_tracts);  
   
  if ~isempty(include)
    mapped_include = find(ismember(bundle_indices,include));
    tracts = tracts(mapped_include,:);
    
    if ignore_bundles
      bundle_indices = include;
    else
      bundle_indices = mapped_include;
    end
    prop_values = prop_values(mapped_include,:);
  end

  %   base_widths = get_properties(prop_keys, prop_values, 'base_width', 1.0, num_tracts);  
  
  num_tracts = size(tracts,1);
  
  max_bundle_index = max(bundle_indices);
  
  if highlight_axes
    num_required_colours = max_bundle_index * 4;  
  else
    num_required_colours = max_bundle_index;
  end
  
  
  set_bundle_colours(num_required_colours);

  %Set up the figure
  if ~hold_on
    main_fig = my_figure(strrep(tracts_filename,'_',''), 1, 1, [1 1 1], 1, [],[],~invisible);
  else
    hold on;
  end  
  cameratoolbar('Show');
  cameratoolbar('SetMode','orbit');
  
  
  % If '-num_width_sections' parameter is not set plot the boundary of the tract
  % instead.
  if strfind(style, 'tracts') ~= 0
    
    if oblong
      error('''-oblong'' option cannot be used with ''tracts'' style');
    end
    
    
    if highlight_axes
      error('''-highlight_axes'' option cannot be used with ''tracts'' style');
    end
    
    add_tracts_to_plot(tracts, colours_of_bundles, ones(num_tracts,1), ones(num_tracts,1), tube_corners, num_length_sections, transparency, bundle_indices);
     
  end
  if strfind(style, 'tubes') ~= 0
  
    [strands, bundle_indices] = tracts2strands(tracts, ones(num_tracts,1), num_width_sections, highlight_axes, oblong, bundle_indices, strands_per_acs, acs);
    tcks = strands2tcks(strands);
    
    radii = ones(size(tcks)) * strand_radius;
    
    add_tcks_to_plot(tcks, radii, colours_of_bundles, bundle_indices); 

    %[tcks, colours_of_bundles] = display_strands(tcks, colours_of_bundles, include, 'bundle', line_style);
    
  end
  if strfind(style, 'lines') ~= 0
    
    [strands, bundle_indices] = tracts2strands(tracts, ones(num_tracts,1), num_width_sections, highlight_axes, oblong, bundle_indices, strands_per_acs, acs);  
    tcks = strands2tcks(strands);    
    
    add_lines_to_plot(tcks, colours_of_bundles, bundle_indices);     
  end  
    
  add_sphere_to_plot(sphere_radius);
        
        
  if any([~isempty(slice_x) ~isempty(slice_y) ~isempty(slice_z)])

      add_mri_slice_to_plot(obs_image, slice_x, slice_y, slice_z)

  elseif ~isempty(num_voxels)

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
    
  daspect ([ 1 1 1 ]);

  light           
  lighting gouraud;


  fprintf('\n');
  disp(['Plotted ' num2str(num_tracts) ' tracts.']);   
  
  
  if colour_key
    add_colour_key(bundle_indices);
  end
  
end
