function A = create_movie_of_samples(directory_name, varargin)


  description = 'Creates a movie from a sample directory';
  
  arguments = {'directory_name', 'The directory name which contains the samples'};

  options = {...
             ...%             'strand_radius  ', 0.02,        'float',  'Radii of the plotted strands';...
            };         
                                    
  supplied_options = parse_options(options, varargin);

  for option_i = 1:size(supplied_options,1)
    if strcmp(supplied_options{option_i,3}, 'string')
      eval([supplied_options{option_i,1} ' = ''' supplied_options{option_i,2} ''';']);
    else
      
      if ~isstr(supplied_options{option_i,2})
        supplied_options{option_i,2} = mat2str(supplied_options{option_i,2});
      end
      
      eval([supplied_options{option_i,1} ' = ' supplied_options{option_i,2} ';']);
    end
  end


  if help_display
    display_help_message(description, arguments, options);
    return;
  end
  
  filenames = list_filenames(directory_name, '0000', '.frr');

  
  
   plot_strands(cell2mat(filenames(1)));
  
   
% 2) Resize the figure window to size of movie required.
set(gcf,'Units','pixels') 
winsize = get(gcf,'Position')

% 
% 3) Record the size of the plot window:
% 

% 
% 4) Adjust size of this window to include the whole figure window (if you require the axes, title and axis labels
%      in the movie):
% 
% >> winsize(1:2) = [0 0];
% 
% 5) Set the number of frames:
% 

 num_samples=length(filenames)

% 
% 6) Create the MATLAB movie matrix:
% 
  A=moviein(num_samples,gcf,winsize)
% 
% 7) Fix the features of the plot window (ensures each frame of the movie is the same size):
% 
%   set(fig1,'NextPlot','replacechildren')
  
  
% 
%  8) Within a loop, plot each picture and save to MATLAB movie matrix:
% 
 for sample_i=1:num_samples
  
   plot_strands(cell2mat(filenames(sample_i)));
   
  % add axis label, legends, titles, etc. in here
   A(:,sample_i)=getframe(gcf,winsize);
   
 end
% 
% This procedure creates a movie stored in a special format that is only readable in MATLAB. The first thing you will want to do is to play the movie:
% 
 movie(gcf,A,30,3,winsize)
% 
% where fig1 is the figure handle, A is the movie matrix, 30 is the number of times to repeat movie, 3 is the number of frames per second and winsize is the size of the movie window. You can also save this movie to a file to be loaded another time or on another machine running MATLAB:
% 
 save([directory_name '.mat'],'A');
% 
% and to reload:
% 
% >> load filename.mat 
%   
%   
%   
%   



end