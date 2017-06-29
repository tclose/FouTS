function plot_peaks(varargin)
%  
% PURPOSE: Plots spikes in TestSpikey probability landscape.
%  
% ARGUMENTS: 
%  
%           spike_filename     Filename containing the spike landscape to plot
%  
% OPTIONS (name, description, type, default):
%  
                
  description = 'Plots peaks in TestSpikey probability landscape.';
  
  arguments = {'peaks_filename', 'Filename containing the spike landscape to plot'};
  
  options = {...
              'image_plot    ', 0,           'bool',     'Image plot instead of surface.';...   
              'invert        ', 0,           'bool',     'Invert the plot upside down.';...  
              'extent        ', 1.5,         'float',    'Extent of the peaks to plot.';... 
              'resolution    ', 0.01,        'float',    'Resolution of the grid';...               
              'transparency  ', 0.5,         'float',    'Transparency of the plot.';...                                         
            };         

  num_args = nargin;
  parse_arguments      
  if (help_display) 
    return;
  end

%   End arguments %

  peaks = load(peaks_filename);
  
  num_peaks = size(peaks,1);
  
  grid_points = -extent:resolution:extent;
  
  Z = zeros(size(X));
  
  for x = grid_points
  
    for y = grid_points
     
      for peak_i = 1:num_peaks
        
      end
      
    end
    
  end
    
  my_figure(peaks_filename);
  
  [X,Y] = meshgrid(grid_points, grid_points);
  
  if image_plot
    
    imagesc(spikes);
    
  else
  
    X = repmat([1:size(spikes,2)]',1, size(spikes,2));
    Y = repmat([1:size(spikes,1)], size(spikes,1), 1);

    h = surf(X,Y,spikes);

    set(h,'FaceAlpha', transparency);
    
    cameratoolbar('Show');
    cameratoolbar('SetMode','orbit');

  end
  
  daspect([1 1 1]);
  
end