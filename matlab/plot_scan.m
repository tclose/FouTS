function plot_scan(varargin)
%  
% PURPOSE: Plots the output of 'scan_strands'
%  
% ARGUMENTS: 
%  
%           scan_filename     The filename of the '.mif', '.sst', or '.tst' scan_image containing the scan.
%  
% OPTIONS (name, description, type, default):
%  
%          -scale_range    
%                 The range of the colour scale plot.
%                 matrix_1x2
%  
%          -state_filename 
%                 Filename of the states to overlay on the surface plot.
%                 string
%                 ''
%  
%          -include_states 
%                 States to include in plot.
%                 matrix_:x:
%  
%          -mark_rand_samples
%                 Mark samples where the momentum is randomized.
%                 bool
%                 0
%  
%          -highlight_samples
%                 Highlight selected samples with a larger marker
%                 matrix_:x:
%  
%          -image          
%                 Plot as a flat image instead of a surface.
%                 bool
%                 0
%  
%          -exponential    
%                 Change plot to display the probability instead of the log probability.
%                 bool
%                 0
%  
%          -offset         
%                 Offset the values by a constant amount (used in conjunction with '-exponential').
%                 float
%                 0
%  
%          -transparency   
%                 The transparency of the surface.
%                 float
%                 0.5

  
  description = 'Plots the output of ''scan_strands''';
  
  arguments = {'scan_filenames', 'The filename of the ''.mif'', ''.sst'', or ''.tst'' scan_image containing the scan.'};

  options = {...
            'scale_range    ', [],          'matrix_1x2', 'The range of the colour scale plot.';...            
            'state_filename ', '',          'string'    , 'Filename of the states to overlay on the surface plot.';...
            'include_states ', [],          'matrix_:x:', 'States to include in plot.';...
            'mark_rand_samples', 1,         'bool',       'Mark samples where the momentum is randomized.';...
            'highlight_samples', [],        'matrix_:x:', 'Highlight selected samples with a larger marker';...
            'image_plot     ', 0,           'bool'      , 'Plot as a flat image instead of a surface.';...
            'map_to_surface ', 0,           'bool'      , 'Maps the probabilities on the surfaces of the scans, representing the probability due to that dimension.';...
            'exponential    ', 0,           'bool'      , 'Change plot to display the probability instead of the log probability.';...
            'offset         ', 0,           'float'     , 'Offset the values by a constant amount (used in conjunction with ''-exponential'').';...            
            'transparency   ', 0.5,         'float'     , 'The transparency of the surface.';...
            };         


  multiple_last_arguments = 1;
          
  num_args = nargin;
  parse_arguments      
  if (help_display) 
    return;
  end

%   End arguments %


  %Expand any directories
  expanded_scan_filenames = cell(0);
  
  for scan_i = 1:length(scan_filenames);
  
    if isdir(scan_filenames{scan_i})
    
      filenames = dir(scan_filenames{scan_i});

      for file_i = 1:length(filenames)
        if strcmp(file_extension(filenames(file_i).name), 'mif')
          expanded_scan_filenames{end+1} = [scan_filenames{scan_i}, filesep, filenames(file_i).name];
        end
        
      end
      
    else      
      expanded_scan_filenames{end+1} = scan_filenames{scan_i};
    end
  
  end
    
  num_scans = length(expanded_scan_filenames);
  

  %If states provided preload them to save having to read the file multiple
  %times.
  if ~isempty(state_filename)

    if strcmp(file_extension(state_filename), 'sst')

      [sets, props, set_prop_keys, set_prop_values]  = load_strand_sets(state_filename,include_states,1);
  
    elseif strcmp(file_extension(state_filename), 'tst') 

      [sets, props, set_prop_keys, set_prop_values]  = load_tract_sets(state_filename,include_states,1);
      
    else
      error(['Unrecognised extension ''' file_extension(state_filename) '''.']);
    end
    
  end
      
  %Plot the scans.
  for scan_i = 1:num_scans
    
    scan_filename = expanded_scan_filenames{scan_i};

    figure_name = strrep(scan_filename,'_',' ');

    if exponential
      figure_name = [figure_name ' - Probability'];
    else
      figure_name = [figure_name ' - Log Probability'];
    end


    %If 1D scan
    if strcmp(file_extension(scan_filename),'sst') || strcmp(file_extension(scan_filename),'tst') || strcmp(file_extension(scan_filename),'str') || strcmp(file_extension(scan_filename),'tct')

      [prop_keys, prop_values] = read_element_properties([scan_filename 'x']);

      scan_data = get_properties(prop_keys, prop_values, 'log_px');

      if offset > 0
        scan_data = scan_data + offset;
      end

      if exponential
        scan_data = exp(scan_data);
      end  

      my_figure(figure_name, scan_i, num_scans);
      
      plot(scan_data);


    %If 2d scan  
    elseif strcmp(file_extension(scan_filename), 'mif')

      scan_image = read_image(scan_filename);

      if ~isfield(scan_image, 'data')
         error(['Could not load scan_image from ' scan_filename ]);
      end

      scan_range1 = [-1:2/(scan_image.dim(1)-1):1];
      scan_range2 = [-1:2/(scan_image.dim(2)-1):1];    

      if offset > 0
        scan_image.data = scan_image.data + offset;
      end

      if exponential
        scan_image.data = exp(scan_image.data);
      end


      %Plot as a flat image.
      if image_plot
        
        my_figure(figure_name, scan_i, num_scans);
        
        if isempty(scale_range)
          imagesc(scan_range1, scan_range2, scan_image.data);
        else
          imagesc(scan_range1, scan_range2, scan_image.data, scale_range);
        end

        colorbar;

      %Else plot as a surface.
      else

        [X,Y] = meshgrid(scan_range1, scan_range2);

        my_figure(figure_name, scan_i, num_scans, [], 1);

        h = surf(X, Y, scan_image.data');


        if transparency ~= 1
          set(h,'FaceAlpha', transparency);
        end

        %If states provided (i.e. a MCMC sequence) plot projection on to
        %scan slice.
        if ~isempty(state_filename)

          axis1 = eval(scan_image.axis1);
          axis2 = eval(scan_image.axis2);
          origin = eval(scan_image.origin);

          if strcmp(file_extension(state_filename), 'sst')

            num_sets = size(sets,1);

            px = get_properties(set_prop_keys, set_prop_values, 'log_px', zeros(num_sets,1), num_sets);

            proj1 = zeros(num_sets,1);
            proj2 = zeros(num_sets,1);  


            for set_i = 1:num_sets

              strands = sets{set_i};

              num_strands = size(strands,1);

              for strand_i = 1:num_strands

                strand_disp = strands{strand_i} - origin{strand_i + 2}{3};

                proj1(set_i) = proj1(set_i) + dot(axis1{strand_i + 2}{3}(:), strand_disp(:));
                proj2(set_i) = proj2(set_i) + dot(axis2{strand_i + 2}{3}(:), strand_disp(:));  

              end

            end


            axis1_norm2 = 0;
            axis2_norm2 = 0;

            for strand_i = 1:num_strands
              axis1_norm2 = axis1_norm2 + axis1{strand_i + 2}{3}(:)' * axis1{strand_i + 2}{3}(:);
              axis2_norm2 = axis2_norm2 + axis2{strand_i + 2}{3}(:)' * axis2{strand_i + 2}{3}(:);
            end

            proj1 = proj1 ./ axis1_norm2;
            proj2 = proj2 ./ axis2_norm2;


          elseif strcmp(file_extension(state_filename), 'tst') 

            num_sets = size(sets,1);

            px = get_properties(set_prop_keys, set_prop_values, 'px', zeros(num_sets,1), num_sets);

            proj1 = zeros(num_sets,1);
            proj2 = zeros(num_sets,1);  


            for set_i = 1:num_sets

              tracts = sets{set_i};

              num_tracts = size(tracts,1);

              for tract_i = 1:num_tracts

                for strand_i = 1:3

                  strand_disp = tracts{tract_i}{strand_i} - origin{tract_i + 2}{strand_i}{3};

                  proj1(set_i) = proj1(set_i) + dot(axis1{tract_i + 2}{strand_i}{3}(:), strand_disp(:));
                  proj2(set_i) = proj2(set_i) + dot(axis2{tract_i + 2}{strand_i}{3}(:), strand_disp(:));  

                end

              end

            end


            axis1_norm2 = 0;
            axis2_norm2 = 0;

            for tract_i = 1:num_tracts

              for strand_i = 1:3
                axis1_norm2 = axis1_norm2 + axis1{tract_i + 2}{strand_i}{3}(:)' * axis1{tract_i + 2}{strand_i}{3}(:);
                axis2_norm2 = axis2_norm2 + axis2{tract_i + 2}{strand_i}{3}(:)' * axis2{tract_i + 2}{strand_i}{3}(:);
              end

            end

            proj1 = proj1 ./ axis1_norm2;
            proj2 = proj2 ./ axis2_norm2;



          else
            error(['Unrecognised extension ''' file_extension(state_filename) '''.']);
          end


          hold on;

          %If map_to_surface option is selected, map probability values via 
          %bilinear interpolation to scan surface.
          if map_to_surface
             
            for set_i = 1:num_sets

              
              if (proj1(set_i) < -1) || (proj1(set_i) > 1) || (proj2(set_i) < -1) || (proj2(set_i) > 1)
              
                px(set_i) = 0;
                
              else
                
                
                index1 = (scan_image.dim(1)-1) * (proj1(set_i) + 1)/2 + 1;                
                index2 = (scan_image.dim(2)-1) * (proj2(set_i) + 1)/2 + 1;                
              
                low1 = floor(index1);
                high1 = ceil(index1);
                frac1 = rem(index1,1);
                
                low2 = floor(index2);
                high2 = ceil(index2);
                frac2 = rem(index2,1);
                
                px(set_i) = (1-frac2) * ((1-frac1) * scan_image.data(low1, low2) + frac1 * scan_image.data(high1, low2)) + ... 
                            frac2 * ((1-frac1) * scan_image.data(low1, high2) + frac1 * scan_image.data(high1, high2));
                
              end

            end
            
            state_plot_colour = [1 0.5 1];
            
          else
            
            state_plot_colour = 'yellow';
            
          end    
          
          
          plot3(proj1, proj2, px, 'Color', state_plot_colour);

          plot3(proj1(1), proj2(1), px(1), 'md', 'MarkerEdgeColor', 'red',  'MarkerFaceColor', 'red', 'MarkerSize', 10);

          if mark_rand_samples

            sample_period = str2double(props.sample_period);

            mark_samples = sample_period:sample_period:num_sets;

            plot3(proj1(mark_samples), proj2(mark_samples), px(mark_samples), 'mo', 'MarkerEdgeColor', 'blue',  'MarkerFaceColor', 'blue', 'MarkerSize', 5);

          end

          if ~isempty(highlight_samples)

            highlight_samples = highlight_samples + 1;

            highlight_samples = highlight_samples(highlight_samples < num_sets);

            plot3(proj1(highlight_samples), proj2(highlight_samples), px(highlight_samples), 'md', 'MarkerEdgeColor', 'green',  'MarkerFaceColor', 'green', 'MarkerSize', 8);

          end

          hold off


        end
        
      end

      axis1_delim = strfind(scan_image.axis1_location,'/');

      if isempty(axis1_delim)
        axis1_basename = 0;
      else
        axis1_basename = scan_image.axis1_location((axis1_delim(end)+1):end);
      end

      xlabel(axis1_basename);


      axis2_delim = strfind(scan_image.axis2_location,'/');

      if isempty(axis2_delim)
        axis2_basename = 0;
      else
        axis2_basename = scan_image.axis2_location((axis2_delim(end)+1):end);
      end

      ylabel(axis2_basename);




    else

      error(['Unrecognised file_extension, ''' file_extension(scan_filename) '''.']);

    end

  end
    
end
