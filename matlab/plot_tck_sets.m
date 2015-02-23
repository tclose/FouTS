function [main_fig, colour_indices] = plot_tcks_sets(varargin)
%  
% PURPOSE: Plots tcks from tck files, and optionally displays reference sphere and voxels
%  
% ARGUMENTS: 
%  
%           tcks_filename     The filename of the tcks in either .tck or .frr formats.
%  
% OPTIONS (name, description, type, default):
%  
%          -tck_radius  
%                 Radii of the plotted tcks
%                 float
%                 0.02
%  
%          -num_points     
%                 Number of points to plot along each tck.
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
%          -tck_include 
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

  description = 'Plots tcks from tck files, and optionally displays reference sphere and voxels';
  
  arguments = {'tcks_filename', 'The filename of the tcks in either .tck or .frr formats.'};
  
  options = {...
            'tck_radius  ', 0.02,        'float',  'Radii of the plotted tcks';...
            'num_points     ', 100,         'int',    'Number of points to plot along each tck.';...            
            'include        ', [],          'matrix_:x:', 'The indices of the sets to include in the plot';...   
            'dont_include        ', [],     'matrix_:x:', 'The indices of the samples not to include in the plot';...  
            'last',            0,           'bool',       'Overrides ''-include'' to only print the last tract set.';...            
            'tck_include ', [],          'matrix_:x:', 'The indices of the tcks to include in the plot';...             
            'bundle_include ', [],          'matrix_:x:', 'The indices of the bundles to include in the plot';...            
            'colours_of_bundles', colours_of_bundles,     'matrix_:x3', 'Colours of the plotted tcks';...
            'voxel_size     ', 0.15,        'float',  'Size of reference voxel';...
            'num_voxels     ', 3,           'int',    'Number of reference voxels';...
            'cube_size      ', 0,           'float',  'Size of reference voxel';...
            'happy_colours',         0,     'bool', 'Uses ''happy colours''';...
            'inv_happy_colours',         0,     'bool', 'Uses the inverse of ''happy colours''';...
            'style          ', 'lines',     'string', 'Plot to style (either ''tubes'' or ''lines'')..';...            
            'tube_corners   ', 6,          'int',    'Change plot to ''line style''.(NB: Only relevant with ''-num_tcks'' option = 0).';...
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
    
  [tck_sets, props, set_prop_keys, set_prop_values, elem_prop_keys, elem_prop_values] = load_tck_sets(tcks_filename, [], 1);
  
  num_sets = size(tck_sets,1);
  
  if last

    if ~isempty(dont_include) || ~isempty(include)
       error('Error! -last cannot be used simultaneously with -include or -dont_include options');
    end
    
    include = size(tck_sets,1)-1;
    
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
  
  if ~isempty(tck_include)
    error('Cannot use ''tck_include'' and ''bundle_include simultaneously');
  end    
  
  include = include + 1;
  
  
  if properties_plot < 2

    overall_max_bundle_index = 0;

    for set_i = 1:num_sets
      
      num_tcks = size(tck_sets{set_i},1);
      
      if size(elem_prop_keys,2)
        prop_values = elem_prop_values{set_i};
      else
        prop_values = [];
      end
      
      bundle_indices = get_properties(elem_prop_keys, prop_values, 'bundle_index', [0:1:(num_tcks-1)]', num_tcks);
      
      max_bundle_index = max(bundle_indices);
      
      if max_bundle_index > overall_max_bundle_index
        overall_max_bundle_index = max_bundle_index;
      end
    end  

    set_bundle_colours(overall_max_bundle_index);    

    %Set up the figure.
    main_fig = my_figure(tcks_filename, 1, 3, [1 1 1], 1, [],[],~invisible);

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

      tcks = tck_sets{set_i};  


      if ~isempty(bundle_include)

        num_tcks = size(tcks,1);

        bundle_indices  = get_properties(elem_prop_keys, prop_values, 'bundle_index', 1:num_tcks);

        tck_include = [];
        
        for tck_i = 1:num_tcks
          if any(find(bundle_include == bundle_indices(tck_i)))
            tck_include = [tck_include; tck_i]; %#ok<AGROW>
          end
        end

        if isempty(tck_include)
          disp('Warning! No bundles matched the given indices, returning all tcks.');
        end

      end

      num_tcks = size(tcks,1);

      bundle_indices  = get_properties(elem_prop_keys, prop_values, 'bundle_index', [0:1:(num_tcks-1)]', num_tcks);
      radii           = get_properties(elem_prop_keys, prop_values, 'tract_radius', tck_radius, num_tcks);


      if ~isempty(tck_include)

        tcks = tcks(tck_include);
        bundle_indices = bundle_indices(tck_include);
        radii = radii(tck_include);

      end
      

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
    disp(['Plotted ' num2str(num_sets) ' tck sets.']);   

  end
  
  if properties_plot > 0 && length(include) ~= 1

    if exist([tcks_filename 'x'])
      
      ext_props_fig = plot_extend_elem_properties([tcks_filename 'x'], set_prop_keys, set_prop_values,include);

    end
  end

  if (properties_plot == 1 || properties_plot == 3) && length(include) ~= 1
    if exist([tcks_filename 'xx'])
      
      ext_ext_props_fig = plot_extended_extend_elem_properties([tcks_filename 'xx'], elem_prop_keys, elem_prop_values, colours_of_bundles, include);
  
    end
    
  end

end