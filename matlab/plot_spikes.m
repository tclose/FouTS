function plot_spikes(varargin)
%  
% PURPOSE: Plots spikes in TestSpikey probability landscape.
%  
% ARGUMENTS: 
%  
%           spike_filename     Filename containing the spike landscape to plot
%  
% OPTIONS (name, description, type, default):
%  
                
  description = 'Plots spikes in TestSpikey probability landscape.';
  
  arguments = {'spikes_filename', 'Filename containing the spike landscape to plot'};
  
  options = {...
              'image_plot    ', 0,           'bool',     'Image plot instead of surface.';...   
              'invert        ', 0,           'bool',     'Inver the plot upside down.';...                           
              'transparency  ', 0.5,         'float',    'Transparency of the plot.';...                                         
            };         

  parse_arguments      
  if (help_display) 
    return;
  end

%   End arguments %

  spikes = load(spikes_filename);
  
  if invert
    spikes = -spikes;
  end

  my_figure(spikes_filename);
  
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