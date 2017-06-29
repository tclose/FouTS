function h = my_figure(name, figure_index, num_figures, data_aspect_ratio, camera_rotate, max_height, max_width, visible)
%function h = my_figure(name, figure_index, num_figures, data_aspect_ratio, camera_rotate, max_height, max_width)
%
%
%Sets up the the figure the way I like it as well as automatically position
%it on the screen.
%h = my_figure(name, figure_index, num_figures, data_aspect_ratio, camera_rotate, max_height, max_width)

  screen_size = get(0,'ScreenSize');
  screen_width = screen_size(3);
  screen_height = screen_size(4);
 
  %Constant references.
  screen_aspect_ratio = screen_width / screen_height;
  title_buffer = 105;
  
    
  %Optional parameters.
  if ~exist('figure_index','var')
    figure_index = 1;
  end

  if ~exist('num_figures','var')
    num_figures = 1;
  end
  
  if ~exist('camera_rotate','var')
    camera_rotate = 0;
  end
  
  if ~exist('data_aspect_ratio','var')
    data_aspect_ratio = [1 1 1];
  end
  
  if ~exist('max_height','var')
    max_height = [];
  end
  
  if isempty(max_height)
      max_height = 0.5 * screen_height;
  end
  
  if ~exist('max_width','var')
    max_width = [];
  end
  
  if isempty(max_width)
      max_width = 0.4 * screen_width;
  end
  
  if ~exist('visible','var')
    visible = 1;
  end
  
 
  if camera_rotate
    title_buffer = title_buffer + 25;
  end
  
  %Validation
  if figure_index > num_figures
    error (['Figure index (' num2str(figure_index) ') exceeds number of figures (' num2str(num_figures) ').']);
  end
  
  fig_position = zeros(1,4);
  
  num_columns = round(sqrt(screen_aspect_ratio * num_figures));
  num_rows = ceil(num_figures / num_columns);

  width = min(max_width , screen_width / num_columns);
  height = min(max_height , screen_height / num_rows - title_buffer);% - title_buffer;

  figure_name = [strrep(name,'_',' ')];

  if visible

    fig_position(3) = width;
    fig_position(4) = height;
  
    if num_columns == 1
      
      fig_position(1) = (screen_width-width)/2;
      
    else
      
      column = mod(figure_index-1, num_columns);
      fig_position(1) = width * column;
      
    end
    
    if num_rows == 1
      
      fig_position(2) = (screen_height-height)/2;
      
    else
      
      row = floor((figure_index-1) / num_columns);
      fig_position(2) = (height + title_buffer) * (num_rows - row - 1);
          
    end  
  
    h = figure();

    set(h,'Units','pixels');
    set(h, 'Position', fig_position);
    set(h, 'Name', figure_name);

  else
    h = figure('Visible','Off');
    set(h, 'Position', [0 0, 500,500]);
  end
  
  whitebg(h,'black'); 

  if ~isempty(data_aspect_ratio)
    daspect(data_aspect_ratio);
  end
  
  if camera_rotate
    set(h, 'Color', 'black');
    set(h, 'DoubleBuffer', 'on');
    cameratoolbar('Show');
    cameratoolbar('SetMode','orbit');
  end
  
  
end
