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
  
  arguments = {'scan_filename', 'The filename of the ''.mif'', ''.sst'', or ''.tst'' scan_image containing the scan.'};

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


  multiple_last_arguments = 0;
          
  parse_arguments      
  if (help_display) 
    return;
  end

%   End arguments %

  figure_name = strrep(scan_filename,'_',' ');

  if exponential
    figure_name = [figure_name ' - Probability'];
  else
    figure_name = [figure_name ' - Log Probability'];
  end


  %If 2d scan  
  if ~strcmp(file_extension(scan_filename), 'mif')
    error(['Unrecognised file_extension, ''' file_extension(scan_filename) '''.']);
  end

  scan_image = read_image(scan_filename);

  if ~isfield(scan_image, 'data')
     error(['Could not load scan_image from ' scan_filename ]);
  end

  axis1_scale = str2num(scan_image.axis1_scale);
  axis2_scale = str2num(scan_image.axis2_scale);
  
  scan_range1 = [-axis1_scale:2*axis1_scale/(scan_image.dim(1)-1):axis1_scale];
  scan_range2 = [-axis2_scale:2*axis2_scale/(scan_image.dim(2)-1):axis2_scale];    

  if offset > 0
    scan_image.data = scan_image.data + offset;
  end

  if exponential
    scan_image.data = exp(scan_image.data);
  end


  %Plot as a flat image.
  if image_plot

    my_figure(figure_name);

    if isempty(scale_range)
      imagesc(scan_range1, scan_range2, scan_image.data);
    else
      imagesc(scan_range1, scan_range2, scan_image.data, scale_range);
    end

    colorbar;

  %Else plot as a surface.
  else

    [X,Y] = meshgrid(scan_range1, scan_range2);

    my_figure(figure_name, 1, 1, [], 1);

    h = surf(X, Y, scan_image.data);


    if transparency ~= 1
      set(h,'FaceAlpha', transparency);
    end

    %If states provided (i.e. a MCMC sequence) plot projection on to
    %scan slice.
    if ~isempty(state_filename)

      trans_data = scan_image.data';
      
      iter_filename = [file_base(state_filename) '.iter.' file_extension(state_filename)];
      
      if exist(iter_filename)
        states = load(iter_filename);
      else
        states = load(state_filename);
      end
      
      props = read_properties([state_filename '.prp']);

      num_states = size(states,1);
      
      proj1 = states(:,1) ./ str2num(scan_image.axis1_scale);
      proj2 = states(:,2) ./ str2num(scan_image.axis2_scale);
      px = zeros(num_states,1);
      
      %If map_to_surface option is selected, map probability values via 
      %bilinear interpolation to scan surface.
      if map_to_surface

        for state_i = 1:num_states

          if (proj1(state_i) >= -1) && (proj1(state_i) <= 1) && (proj2(state_i) >= -1) && (proj2(state_i) <= 1)

            index1 = (scan_image.dim(1)-1) * (proj1(state_i) + 1)/2 + 1;                
            index2 = (scan_image.dim(2)-1) * (proj2(state_i) + 1)/2 + 1;                

            low1 = floor(index1);
            high1 = ceil(index1);
            frac1 = rem(index1,1);

            low2 = floor(index2);
            high2 = ceil(index2);
            frac2 = rem(index2,1);

            px(state_i) = (1-frac2) * ((1-frac1) * trans_data(low1, low2) + frac1 * trans_data(high1, low2)) + ... 
                        frac2 * ((1-frac1) * trans_data(low1, high2) + frac1 * trans_data(high1, high2));

          end

        end

        state_plot_colour = [1 0.5 1];

      else

        state_plot_colour = 'yellow';

      end    
      

      hold on;


      plot3(states(:,1), states(:,2), px, 'Color', state_plot_colour);

      plot3(states(1,1), states(1,2), px(1), 'md', 'MarkerEdgeColor', 'red',  'MarkerFaceColor', 'red', 'MarkerSize', 10);

      if mark_rand_samples

        sample_period = str2double(props.sample_period);

        mark_samples = sample_period:sample_period:num_states;

        plot3(states(mark_samples,1), states(mark_samples,2), px(mark_samples), 'mo', 'MarkerEdgeColor', 'blue',  'MarkerFaceColor', 'blue', 'MarkerSize', 5);

      end

      if ~isempty(highlight_samples)

        highlight_samples = highlight_samples + 1;

        highlight_samples = highlight_samples(highlight_samples < num_states);

        plot3(states(highlight_samples,1), states(highlight_samples,2), px(highlight_samples), 'md', 'MarkerEdgeColor', 'green',  'MarkerFaceColor', 'green', 'MarkerSize', 8);

      end

      hold off

    end

  end

  xlabel('Axis 1');
  ylabel('Axis 2');

    
end
