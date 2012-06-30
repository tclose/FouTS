function plot_explore_stats(varargin)
%  
% PURPOSE: Plots the output statistics of a parameter-space exploration
%  
% ARGUMENTS: 
%  
%           dir_name     The name of the directory that the samples are stored in.
%  
% OPTIONS (name, description, type, default):
%

  description = 'Plots the output statistics of a parameter-space exploration';
  
  arguments = {'dir_name', 'The name of the directory that the samples are stored in.'};

  options = {...
                 'chunk_to_plot',   1,           'int',    'Index of the chunk which will be plotted';...
...%             'strand_radius  ', 0.02,        'float',  'Radii of the plotted strands';...
            };         

  % Call a script called 'parse_arguments' to parse all arguments and
  % options and set the desired variables.
  
  parse_arguments      
  if (help_display) return; end
         
  stats_types = {...
    'strong_dist';...
    'strong_px';...
    ...%'strong_central_px';...
    'weak_dist';...
    'weak_px';...
    ...%'weak_central_px';...
    ...%'max_arc_top_dist'...
  }  

  stats = [];

  for stats_i = 1:length(stats_types)
     
    loaded_stats = load([dir_name '/' stats_types{stats_i} '.txt']);
    
    if size(loaded_stats,1) > 10000    
      loaded_stats = loaded_stats(1:10000);
    end
    
    stats = [stats reshape(loaded_stats, [size(loaded_stats,1) 1 size(loaded_stats,2)])];
    
  end
  
  reshape_vector = [size(stats,1) * size(stats,3) , 1];
  
  f = figure;
  set(f, 'Name', 'Strong');
  
  scatter(reshape(stats(:,1,:),reshape_vector), reshape(stats(:,2,:), reshape_vector));
  whitebg('black');

  g = figure;
  set(g, 'Name', 'Weak');
  scatter(reshape(stats(:,3,:),reshape_vector), reshape(stats(:,4,:), reshape_vector));
  whitebg('black');  
  
  
  figure;
  plot(stats(:,:,chunk_to_plot));
  whitebg('black');
  legend(strrep(stats_types(:), '_', ' '));
%     legend('distance', 'full-strong', 'outside-weak', 'outside-strong','central-voxel', 'b0');
  legend('Location', 'EastOutside');
end
