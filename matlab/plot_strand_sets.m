function main_fig = plot_strands_sets(varargin)
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
%                 The indices of the sets to include in the plot
%                 matrix_:x:
%  
%          -last
%                 Overrides '-include' to only print the last tract set.
%                 bool
%                 0
%  
%          -strand_include 
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
%                 Plot to style (either 'tubes' or 'lines')..
%                 string
%                 'lines'
%  
%          -properties_plot
%                 Do, don't or only plot the extended properties associated with the tract set (0 -> don't plot properties, 1 -> plot properties, 2 -> only plot properties).
%                 int
%                 1
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

  main_fig = -1;

  global colours_of_bundles;

  description = 'Plots strands from strand files, and optionally displays reference sphere and voxels';
  
  arguments = {'strands_filename', 'The filename of the strands in either .tck or .frr formats.'};
  
  options = {...
            'strand_radius  ', 0.02,        'float',  'Radii of the plotted strands';...
            'num_points     ', 100,         'int',    'Number of points to plot along each strand.';...            
            'include        ', [],          'matrix_:x:', 'The indices of the sets to include in the plot';...   
            'dont_include        ', [],     'matrix_:x:', 'The indices of the samples not to include in the plot';...  
            'last',            0,           'bool',       'Overrides ''-include'' to only print the last tract set.';...            
            'strand_include ', [],          'matrix_:x:', 'The indices of the strands to include in the plot';...             
            'bundle_include ', [],          'matrix_:x:', 'The indices of the bundles to include in the plot';...            
            'colours_of_bundles', colours_of_bundles,     'matrix_:x3', 'Colours of the plotted strands';...
            'voxel_size     ', 0.15,        'float',  'Size of reference voxel';...
            'num_voxels     ', 3,           'int',    'Number of reference voxels';...
            'cube_size      ', 0,           'float',  'Size of reference voxel';...
            'happy_colours',         0,     'bool', 'Uses ''happy colours''';...
            'inv_happy_colours',         0,     'bool', 'Uses the inverse of ''happy colours''';...
            'style          ', 'lines',     'string', 'Plot to style (either ''tubes'' or ''lines'')..';...            
            'tube_corners   ', 6,          'int',    'Change plot to ''line style''.(NB: Only relevant with ''-num_strands'' option = 0).';...
            'properties_plot', 0,           'int',    'Do, don''t or only plot the extended properties associated with the tract set (0 -> don''t plot properties, 1 -> plot properties, 2 -> only plot properties).';...
            'sphere_radius  ', 0,           'float',  'Size of reference sphere';...
            'no_axes',       0,             'bool', 'Removes axes from plot';...
            'clean',         0,             'bool', 'Removes everything else from plot';...
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
    
  [strand_sets, props, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values] = load_strand_sets(strands_filename, [], 1);
  
  num_sets = size(strand_sets,1);
  
  if last

    if ~isempty(dont_include) || ~isempty(include)
       error('Error! -last cannot be used simultaneously with -include or -dont_include options');
    end
    
    include = size(strand_sets,1)-1;
    
  else
   
    if ~isempty(dont_include)

      if ~isempty(include)
        error('Error! -include and -dont_include options cannot be used simultaneously');
      end

      include = 1:num_sets;
      
      include(dont_include) = [];
      
      include = include-1;

    end
  
  end
  
  if ~isempty(strand_include)
    error('Cannot use ''strand_include'' and ''bundle_include simultaneously');
  end    
  
  include = include + 1;
  
  
  if properties_plot < 2

    overall_max_bundle_index = 0;

    for set_i = 1:num_sets
      
      num_strands = size(strand_sets{set_i},1);
      
      if size(elem_prop_keys,2)
        prop_values = elem_prop_values{set_i};
      else
        prop_values = [];
      end
      
      bundle_indices = get_properties(elem_prop_keys, prop_values, 'bundle_index', [0:1:(num_strands-1)]', num_strands);
      
      max_bundle_index = max(bundle_indices);
      
      if max_bundle_index > overall_max_bundle_index
        overall_max_bundle_index = max_bundle_index;
      end
    end  

    set_bundle_colours(overall_max_bundle_index);    

    %Set up the figure.
    main_fig = my_figure(strands_filename, 1, 3, [1 1 1], 1, [],[],~invisible);

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
    
    daspect ([ 1 1 1 ]);
   
    add_sphere_to_plot(sphere_radius, props);

    if voxel_size ~= 0

      add_vox_lines_to_plot(voxel_size,num_voxels,~no_voxline_highlight);

    end


    light           
    lighting gouraud;  

    
    if isempty(include)
      include = 1:num_sets;
    end

    for set_i = include
      
      if size(elem_prop_keys,2)
        prop_values = elem_prop_values{set_i};
      else
        prop_values = [];
      end

      strands = strand_sets{set_i};  


      if ~isempty(bundle_include)

        num_strands = size(strands,1);

        bundle_indices  = get_properties(elem_prop_keys, prop_values, 'bundle_index', 1:num_strands);

        strand_include = [];
        
        for strand_i = 1:num_strands
          if any(find(bundle_include == bundle_indices(strand_i)))
            strand_include = [strand_include; strand_i]; %#ok<AGROW>
          end
        end

        if isempty(strand_include)
          disp('Warning! No bundles matched the given indices, returning all strands.');
        end

      end

      num_strands = size(strands,1);

      bundle_indices  = get_properties(elem_prop_keys, prop_values, 'bundle_index', [0:1:(num_strands-1)]', num_strands);
      radii           = get_properties(elem_prop_keys, prop_values, 'tract_radius', strand_radius, num_strands);


      if ~isempty(strand_include)

        strands = strands(strand_include);
        bundle_indices = bundle_indices(strand_include);
        radii = radii(strand_include);

      end
      
      
      tcks = strands2tcks(strands);

      if strfind('tubes', style) == 1

        add_tcks_to_plot(tcks, radii, colours_of_bundles, bundle_indices, tube_corners); 

      elseif strfind('lines', style) == 1

        add_lines_to_plot(tcks, colours_of_bundles, bundle_indices);    

      else

        error(['Urecogised style option ''' style '''.']);

      end

      daspect ([ 1 1 1 ]);

    end

    fprintf('\n');
    disp(['Plotted ' num2str(num_sets) ' strand sets.']);   

  end
  
  if properties_plot > 0 && length(include) ~= 1

    if exist([strands_filename 'x'])
      
      ext_props_fig = plot_extend_elem_properties([strands_filename 'x'], set_prop_keys, set_prop_values,include);

    end
  end

  if (properties_plot == 1 || properties_plot == 3) && length(include) ~= 1
    if exist([strands_filename 'xx'])
      
      ext_ext_props_fig = plot_extended_extend_elem_properties([strands_filename 'xx'], elem_prop_keys, elem_prop_values, colours_of_bundles, include);
  
    end
    
  end

end