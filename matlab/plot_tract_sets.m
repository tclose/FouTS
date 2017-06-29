function [main_fig, num_voxels, voxel_lengths, voxel_offsets, colour_indices] = plot_tracts_sets(varargin)
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
%                 The indices of the samples to include in the plot
%                 matrix_:x:
%  
%          -last
%                 Overrides '-include' to only print the last tract set.
%                 bool
%                 0
%  
%          -tract_include  
%                 The indices of the tracts to include in the plot
%                 matrix_:x:
%  
%          -colours
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
%          -strand_radius  
%                 Radii of the plotted strands. (NB: can only be used with '-style' option 'strands').
%                 float
%                 0.02
%  
%          -style          
%                 Change plot to 'line style'.(NB: Only relevant with '-num_width_sections' option > 0).
%                 string
%                 'tracts'
%  
%          -tube_corners   
%                 Set the number of tube corners for 'tracts' style.
%                 int
%                 12
%  
%          -num_length_sections
%                 Number of samples along the length of the tracts.
%                 int
%                 100
%  
%          -properties_plot
%                 Do, don't or only plot the extended properties associated with the tract set (0 -> don't plot properties, 1 -> plot properties, 2 -> only plot properties).
%                 int
%                 1
%  
%          -new_plot
%                 Plots strands in new figures for comparision
%                 bool
%                 0
%  
%          -transparency
%                 Fix the transparency of the tracts to a set value.  If 0 the intensities of the of the tracts are used instead.
%                 float
%                 1
%  
%          -sphere_radius  
%                 Size of reference sphere
%                 float
%                 0

  main_fig = -1;

  description = 'Plots strands from strand files, and optionally displays reference sphere and voxels';
  
  arguments = {'tracts_filename', 'The filename of the tracts (''.tct'' format).'};
  
  options = {...
            'include        ', [],          'matrix_:x:', 'The indices of the samples to include in the plot';...  
            'dont_include        ', [],     'matrix_:x:', 'The indices of the samples not to include in the plot';...  
            'last',            0,           'bool',       'Overrides ''-include'' to only print the last tract set.';...
            'tract_include  ', [],          'matrix_:x:', 'The indices of the tracts to include in the plot';... 
            'colours', [],     'matrix_:x3', 'Colours of the plotted strands';...
            'colour_indices',           [],         'matrix_:x1',   'The mapping from bundle indices to the colours used';... 
            'compact_colours',  0,          'bool',   'Compact the range of colours created to only use the number needed';...
            'num_voxels     ',          [],         'matrix_1x3',   'Number of reference voxels';...
            'voxel_lengths',            [],         'matrix_1x3',   'Size of reference voxel';...
            'voxel_offsets  ',          [],         'matrix_1x3',   'Offset of reference voxels';...
            'voxel_transparency',       0.25,       'float',        'The alpha value of the voxel lines';...
            'cube_size      ', 0,           'float',  'Size of reference voxel';...
            'strand_radius  ', 0.02,        'float',  'Radii of the plotted strands. (NB: can only be used with ''-style'' option ''strands'').';...
            'style          ', 'tracts',    'string', 'Change plot to ''line style''.(NB: Only relevant with ''-num_width_sections'' option > 0).';...
            'tube_corners   ', 12,          'int',    'Set the number of tube corners for ''tracts'' style.';...
            'num_length_sections', 100,      'int',    'Number of samples along the length of the tracts.';...
            'properties_plot', 0,           'int',    'Do, don''t or only plot the extended properties associated with the tract set (0 -> don''t plot properties, 1 -> plot properties, 2 -> only plot properties).';...
            'new_plot',      0,             'bool',   'Plots strands in new figures for comparision';...
            'transparency'   ,1,            'float',   'Fix the transparency of the tracts to a set value.  If 0 the intensities of the of the tracts are used instead.';...
            'sphere_radius  ', 0,    'float',  'Size of reference sphere';...
            'no_axes',        0,     'bool', 'Removes axes from plot';...
            'num_width_sections    ', -1,           'int',    'Number of strands to plot along each axis. If zero the default 3D surface option will be used instead.';...            
            'strands_per_acs', -1, 'float', 'Instead of a fixed number of width sections the number of plotted strands is determined by the ACS of the tract';...
            'clean',         0,             'bool', 'Removes everything else from plot';...
            'no_voxline_highlight', 0, 'bool', 'Doesn''t highlight corner axes of voxel lines.';...
            'oblong',         0,     'bool', 'Prints an oblong tractlet rather than cylindrical';...
            'happy_colours',         0,     'bool', 'Uses ''happy colours''';...
            'inv_happy_colours',         0,     'bool', 'Uses the inverse of ''happy colours''';...
            'true_tracts_plot',    0,     'bool', 'Plots true tracts (as determined from observed image)';...
            'no_density_plot',    0,     'bool', 'Omitts density from properties plot';...
            'hold_on',        0, 'bool', 'Adds the plot to the previous figure';...
            'highlight_axes', 0,     'bool', 'Highlights the axes when printing in ''strand'' or ''line'' style';...
            'invisible', 0,     'bool', 'Makes the figure invisible (for automatically saving afterwards).';...
            'obs_image_loc', [],   'string', 'Overlays a slice along the x observed image.';...
            'slice_x', [],     'matrix_:x:', 'Overlays slices of observed image along the given indices.';...
            'slice_y', [],     'matrix_:x:', 'Overlays slices of observed image along the given indices.';...
            'slice_z', [],     'matrix_:x:', 'Overlays slices of observed image along the given indices.'};         

  args = parse_arguments(description, arguments, options, varargin);
  if (args == 'help') 
    return;
  end

  if clean
    voxel_size = 0;
    sphere_radius = 0;
    no_axes = 1;
  end
  
  
% End arguments %
  colours = get_happy_colours(colours, happy_colours, inv_happy_colours);
  
  if ~strcmp(file_extension(tracts_filename), 'tst')
    error (['Extension, ''' file_extension(tracts_filename) ''' is not a valid tract set file (''.tst'').']);
  end
  
  
  if (highlight_axes && ~mod(num_width_sections,2))
    error(['If ''-highlight_axes'' option is selected the number of width sections ''-num_width_sections'' needs to be odd (' num2str(num_width_sections) ' provided)'])
  end
  
  [tract_sets, properties, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values] = load_tract_sets(tracts_filename); 
  
  num_sets = size(tract_sets,1);
  
  if isempty(num_voxels)    %#ok<NODEF>
    [num_voxels, true_location, voxel_lengths, voxel_offsets, obs_image] = get_observed_properties(properties, obs_image_loc);
    if ~isempty(num_voxels)
        num_voxels = num_voxels(1:3);
        voxel_lengths = voxel_lengths(1:3);
    end
  elseif isempty(voxel_offsets) && all(num_voxels ~= 0) %#ok<NODEF>
      voxel_offsets = - num_voxels .* voxel_lengths ./2;
  end
  
  
  if (ndims(num_voxels) == 1 || true_tracts_plot)
      if (num_voxels < 0 || true_tracts_plot) 
        num_voxels = num_voxels(1:3);
        voxel_lengths = voxel_lengths(1:3);
      end
  end  

  [loaded_num_length_sections, loaded_num_width_sections] = get_num_section_properties(properties,style);
  
  if num_length_sections == -1
    num_length_sections = loaded_num_length_sections;
  end
  
  if num_width_sections == -1
    num_width_sections = loaded_num_width_sections;
  end
  

  if last

    if ~isempty(dont_include) || ~isempty(include)
       error('Error! -last cannot be used simultaneously with -include or -dont_include options');
    end
    
    include = size(tract_sets,1)-1;
    
  else
   
    if ~isempty(dont_include)

      if ~isempty(include)
        error('Error! -include and -dont_include options cannot be used simultaneously');
      end

      include = [1:num_sets];
      
      include(dont_include) = [];
      
      include = include -1;

    end
  
  end
     
  if length(include) == 1
    properties_plot = 0;
  end
  
  if properties_plot == 0
    num_figs = 1;
  elseif properties_plot == 1
    num_figs = 3;
  elseif properties_plot == 2
    num_figs = 2;
  end
    
  if true_tracts_plot
    num_figs = num_figs + 1;
  end
  
  if properties_plot ~= 2
     
    all_bundle_indices = [];

    for set_i = 1:num_sets
      
      num_tracts = size(tract_sets{set_i},1);
      
      if size(elem_prop_keys,2)
        prop_values = elem_prop_values{set_i};
      else
        prop_values = [];
      end
      
      bundle_indices = get_properties(elem_prop_keys, prop_values, 'bundle_index', [0:1:(num_tracts-1)]', num_tracts);
      
      all_bundle_indices = sort(unique([all_bundle_indices; bundle_indices]));
    end  

    [colours, colour_indices] = set_bundle_colours(all_bundle_indices, colours, colour_indices, compact_colours, highlight_axes); %#ok<NODEF>

    %Set up the figure.
    if ~hold_on
      main_fig = my_figure(tracts_filename, 1, num_figs, [1 1 1], 1, [],[],~invisible);
    else
      main_fig = gcf;
    end
        
    if isempty(include)
      include = 1:num_sets;
    else
      include = include+1; %Convert to Matlab indexing
    end

    for set_i = include

      prop_values = elem_prop_values{set_i};
      tracts = tract_sets{set_i};

      if ~isempty(tract_include)
        tracts = tracts(tract_include+1,:);
        prop_values = prop_values(tract_include+1,:);
      end

      num_tracts = size(tracts,1);  

      intensities = get_properties(elem_prop_keys, prop_values, 'intensity', 1.0, num_tracts);
      acs = get_properties(elem_prop_keys, prop_values, 'alpha', 1.0, num_tracts).^2;
      bundle_indices = get_properties(elem_prop_keys, prop_values, 'bundle_index', (0:1:(num_tracts-1))', num_tracts);

      % If the hold_on option is used, then the tracts are added to the
      % previous plot for reference
      if hold_on
        hold on
      end
      
      num_plots = 0;
      
      if strfind(style, 'tracts') ~= 0

        add_tracts_to_plot(tracts, colours, intensities, ones(num_tracts,1),...
                           tube_corners, num_length_sections, transparency,...
                           bundle_indices, colour_indices);

        num_plots = num_plots + 1;
        
      end
      if strfind(style, 'tubes') ~= 0

        [strands, bundle_indices] = tracts2strands(tracts, ones(num_tracts,1), num_width_sections, highlight_axes, oblong, bundle_indices, strands_per_acs, acs, colour_indices);
        tcks = strands2tcks(strands, num_length_sections);

        radii = ones(size(tcks)) * strand_radius;

        add_tcks_to_plot(tcks, radii, colours, bundle_indices, tube_corners, 1, colour_indices); 

        num_plots = num_plots + 1;
        
      end
      if strfind(style, 'lines') ~= 0

        [strands, bundle_indices] = tracts2strands(tracts, ones(num_tracts,1), num_width_sections, highlight_axes, oblong, bundle_indices, strands_per_acs, acs, colour_indices);  
        tcks = strands2tcks(strands, num_length_sections);    

        add_lines_to_plot(tcks, colours, bundle_indices, '-', colour_indices);     

        num_plots = num_plots + 1;
      end
      
      if num_plots == 0

        error(['Urecogised style option ''' style '''.']);

      end  

    end  

    add_sphere_to_plot(sphere_radius, properties);
        
    if any([~isempty(slice_x) ~isempty(slice_y) ~isempty(slice_z)])
        
        add_mri_slice_to_plot(obs_image, slice_x, slice_y, slice_z)
        
    elseif ~isempty(num_voxels) && all(num_voxels ~= 0)

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
    disp(['Plotted ' num2str(length(include)) ' of ' num2str(num_sets) ' tract sets.']);   
    
  end
  
  % Switch hold off again for the properties plots
  if hold_on
    hold off
  end
  
  if properties_plot && length(include) ~= 1

    if exist([tracts_filename 'x'])

      ext_props_fig = plot_extend_properties([tracts_filename 'x'], set_prop_keys, set_prop_values,include, num_figs - true_tracts_plot - 1, num_figs);

    end


    if exist([tracts_filename 'xx'])
      
      if no_density_plot
%         omit_properties = cell(1,1)
        omit_properties{1} = 'density';
      else
        omit_properties = cell(0);
      end
      
      
      ext_ext_props_fig = plot_extend_elem_properties([tracts_filename 'xx'], elem_prop_keys, elem_prop_values, colours, include, tract_include, num_figs - true_tracts_plot, num_figs, omit_properties);
  
    end
    
  end

  
  if true_tracts_plot
    
    if isempty(true_location)
      error('Could not find true state location in observed image.');
    end
   
    my_figure(['True config: ' true_location], num_figs, num_figs);
    
    cameratoolbar('Show');
    cameratoolbar('SetMode','orbit');
    
    
    if (file_extension(true_location) == 'tck')
    
      [true_tcks, ~, true_prop_keys, true_prop_values] = load_tcks(true_location);
      
      num_true_tcks = size(true_tcks,1);
      
      true_bundle_indices = get_properties(true_prop_keys, true_prop_values, 'bundle_index', [0:1:(num_true_tcks-1)]', num_true_tcks);
      true_radii = get_properties(true_prop_keys, true_prop_values, 'track_radius', 0.03, num_true_tcks);
      
      %Map bundle indices back consecutive range from 0 upwards.
      mapped_indices = true_bundle_indices;
      unique_indices = unique(true_bundle_indices);
      for unique_i = 1:length(unique_indices)
        mapped_indices(true_bundle_indices == unique_indices(unique_i)) = unique_i - 1;
      end
      
      set_bundle_colours(max(unique_indices));    
      
      add_tcks_to_plot(true_tcks,true_radii,colours, mapped_indices,tube_corners);
      
    elseif (file_extension(true_location) == 'tct')
    
      true_tracts = load_tracts(true_location);

      num_true_tracts = size(true_tracts,1);

      set_bundle_colours(num_true_tracts-1);

      add_tracts_to_plot(true_tracts, colours, ones(num_true_tracts,1), ones(num_true_tracts,1), tube_corners, num_length_sections, 1.0);

    else
      error(['Unrecognised extension of true location file ''' true_location '''.']);
    end
    
    set(gcf, 'color', [0 0 0]);

    label_handle = get(gca,'xlabel');
    set(label_handle, 'string', 'X-axis', 'color', [1 1 1]);

    label_handle = get(gca,'ylabel');
    set(label_handle, 'string', 'Y-axis', 'color', [1 1 1]);

    label_handle = get(gca,'zlabel');
    set(label_handle, 'string', 'Z-axis', 'color', [1 1 1]);

    daspect ([ 1 1 1 ]);

    light           
    lighting gouraud;
    
  end
  
  
end