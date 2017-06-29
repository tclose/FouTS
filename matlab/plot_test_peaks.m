function plot_test_peaks(varargin)
%  
% PURPOSE: Plots probability landscape of Prob::TestPeaks function.
%  
% ARGUMENTS: 
%  
%           peaks_filename     Filename containing the spike landscape to plot
%  
% OPTIONS (name, description, type, default):
%  
%          -image_plot
%                 Image plot instead of surface.
%                 bool
%                 0
%  
%          -invert
%                 Invert the plot upside down.
%                 bool
%                 0
%  
%          -extent_of_plot
%                 Extent of the peaks to plot.
%                 float
%                 15
%  
%          -roi_radius
%                 The region of interest radius, after which the pdf drops off with the square of the distance
%                 float
%                 15
%  
%          -barrier_rate
%                 The rate at which it drops off
%                 float
%                 0.3
%  
%          -resolution
%                 Resolution of the grid
%                 float
%                 1
%  
%          -transparency
%                 Transparency of the plot.
%                 float
%                 0.5
%  
%          -samples_location
%                 Location of the samples to plot
%                 string
%                 ''
%  
%          -style_of_samples
%                 Style in which the samples will be plotted.
%                 string
%                 '-'
%  
%          -map_samples
%                 Map samples onto probability landscape.
%                 bool
%                 0
%  
%          -start_sample
%                 Highlight the start sample.
%                 bool
%                 0
%  
%          -mark_every_sample
%                 Mark every nth sample along path (only relevant for iteration plot).
%                 int
%                 0
%  
%          -include_samples
%                 Samples to include in the plots.
%                 matrix_1x:
                
  description = 'Plots probability landscape of Prob::TestPeaks function.';
  
  arguments = {'peaks_filename', 'Filename containing the spike landscape to plot'};
  
  options = {...
              'image_plot'      , 0,           'bool',      'Image plot instead of surface.';...   
              'invert'          , 0,           'bool',      'Invert the plot upside down.';...  
              'extent_of_plot'  , 15,          'float',     'Extent of the peaks to plot.';... 
              'roi_radius'      , 15,          'float',     'The region of interest radius, after which the pdf drops off with the square of the distance';...
              'barrier_rate'    , 0.3,         'float',     'The rate at which it drops off';...              
              'resolution'      , 1,           'float',     'Resolution of the grid';...               
              'transparency'    , 0.5,         'float',     'Transparency of the plot.';...                                         
              'samples_location', '',          'string',    'Location of the samples to plot';...
              'style_of_samples', '-',         'string',    'Style in which the samples will be plotted.';...   
              'map_samples'     , 0,           'bool',      'Map samples onto probability landscape.';...
              'start_sample'    , 0,           'bool',      'Highlight the start sample.';...              
              'mark_every_sample', 0,          'int',       'Mark every nth sample along path (only relevant for iteration plot).';...                            
              'include_samples' , [],          'matrix_1x:','Samples to include in the plots.';...                                
            };         

  num_args = nargin;
  parse_arguments      
  if (help_display) 
    return;
  end
  
%   End arguments %

  peaks = load(peaks_filename);
  
  grid_points = [-extent_of_plot:resolution:extent_of_plot]';
  
  num_grid_points = size(grid_points,1);
  
  Prob = zeros(num_grid_points, num_grid_points);
  
  [X,Y] = meshgrid(grid_points, grid_points);
  
  for x = 1:num_grid_points

    for y = 1:num_grid_points

      point = [grid_points(x) grid_points(y)];

      Prob(x,y) = test_peak_log_prob(point, peaks, roi_radius, barrier_rate);

    end

  end
    
  my_figure(peaks_filename);
  

  
  h = surf(X,Y,Prob);
  
  set(h,'FaceAlpha', transparency);

  
  cameratoolbar('Show');
  cameratoolbar('SetMode','orbit');
  
  
  
  if ~isempty(samples_location)
    
    samples = load(samples_location);

    num_samples = size(samples,1);    
    
    sample_probs = zeros(num_samples,1);
    
    if isempty(include_samples)
      include_samples = 1:num_samples;
    end
    
    if map_samples

      sample_probs = test_peak_log_prob(samples(include_samples,:), peaks, roi_radius, barrier_rate);
      
    end
      
    hold on;
    
    plot3(samples(include_samples,2), samples(include_samples,1), sample_probs,style_of_samples);% 'mo', 'MarkerEdgeColor', 'blue',  'MarkerFaceColor', 'blue', 'MarkerSize', 4);
    
    hold off;
    
    if start_sample
      hold on;
  
      if map_samples      
        start_prob = test_peak_log_prob(samples(1,:), peaks, roi_radius, barrier_rate);
      else
        start_prob = 0;
      end

      plot3(samples(1,2), samples(1,1), start_prob, 'md', 'MarkerEdgeColor', 'red',  'MarkerFaceColor', 'red', 'MarkerSize', 10);
      
        
      hold off;
    end
    
    if mark_every_sample
      
      if map_samples      
        mark_probs = test_peak_log_prob(samples(1:mark_every_sample:num_samples,:), peaks, roi_radius, barrier_rate);
      else
        mark_probs = zeros(size([1:mark_every_sample:num_samples]'));
      end
      
      hold on;
      
      plot3(samples(1:mark_every_sample:num_samples,2), samples(1:mark_every_sample:num_samples,1), mark_probs, 'mo', 'MarkerEdgeColor', 'blue',  'MarkerFaceColor', 'blue', 'MarkerSize', 4);      
      
      hold off;
      
    end
    
  end  
  
end