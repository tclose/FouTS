function colours_of_evolved_strands = plot_evolution(dir_name, varargin)
% 
% 
% PURPOSE: Plots samples from an MCMC sampling of multiple strands described by fourier components
%  
% ARGUMENTS: 
%  
%           collated_samples_dir     The directory that contains the collated samples (output from 'collate samples' program).
%  
% OPTIONS (name, description, type, default):
%  
%          -true_strands_path     The path to a .frr file containing the positions of the true strands
%                              string
%                              ''
%  
%          -include_samples       The indices of the samples to include
%                              matrix_:x:
%                              []
%  
%          -include_strands       The indices of the strands to include
%                              matrix_:x:
%                              []
%  
%          -colours_of_evolved_strands               Base colours for samples of each strand
%                              matrix_:x3
%                              [0.228801116693986 0.363610266029849 0.983863430198537;0.616023151232531 0.815277162004599 0.0953236751558441;0.51015054259208 0.168093971458644 0.739270506551875;0.350622850963505 0.716324475389521 0.294423983509783]
%  
%          -start_shade           Starting shade of the sample colours
%                              float
%                              0.2
%  
%          -end_shade             Ending shade of the sample colours
%                              float
%                              0.8
%  
%          -line_style            Show samples as an evolving line instead of plusses.
%                              bool
%                              0
%  
%          -invert_shades         Invert shades on plotted colours so that shades approach white.
%                              bool
%                              0

  global colours_of_evolved_strands;


  file_paths = list_filenames(dir_name, 'degree_', '.tck');

	if length(file_paths) == 0
		error(['No files loaded from directory ' dir_name '.']);
	end

	strands = load_strands(file_paths{1});
	
	% If the number of previously generated strand colours is less than the number of strands, then reset it and generate again.	
	if (size(strands,1) - 1) > size(colours_of_evolved_strands,1)
		colours_of_evolved_strands = [];
	end


  description = 'Plots samples from an MCMC sampling of multiple strands described by fourier components';
  
  arguments = {'collated_samples_dir', 'The directory that contains the collated samples (output from ''collate samples'' program).'};

  options = {...
            'true_strands_path  ',  [],   'string',      'The path to a .frr file containing the positions of the true strands';...
            'include_samples    ', [],    'matrix_:x:',  'The indices of the samples to include';...            
            'include_strands    ', [],    'matrix_:x:',  'The indices of the strands to include';...            
            'colours_of_evolved_strands            ', colours_of_evolved_strands,    'matrix_:x3',  'Base colours for samples of each strand';...
            'start_shade        ', 0.2,  'float',       'Starting shade of the sample colours';...
            'end_shade          ', 1.0,   'float',       'Ending shade of the sample colours';...
            'line_style         ', 0,     'bool',        'Show samples as an evolving line instead of plusses.';...            
            'invert_shades      ', 0,     'bool',        'Invert shades on plotted colours so that shades approach white.'...            
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

  
  if ~isempty(true_strands_path)
    true_strands = load_strands(true_strands_path);
  else
    true_strands = [];    
  end
  
  if line_style
    style = '-';
    marker = 'none';
  else
    style = 'none';
    marker = '+';
  end
  
  
  shade_range = end_shade - start_shade;


  
  for degree_i = 1:size(file_paths,1) 
    
    strands = load_strands(file_paths{degree_i});
    
    title_name = strrep(strands{size(strands,1),1},'_', ' ');
    
    strands = strands(1:(size(strands,1)-1), :);
    
    display_colours = 0;
    
    if isempty(colours_of_evolved_strands)
      
      display_colours = 1;
      
      for strand_i = 1:size(strands,1)
        colour = rand(1,3);
        
        if invert_shades
          colour = [1 1 1] - colour ./ norm(colour);
        else
          colour = colour ./ norm(colour);
        end
        
        colours_of_evolved_strands = [colours_of_evolved_strands ; colour];
      end

    end
    
    if (isempty(include_strands))
      include_strands = [1:size(strands,1)];
    end
    
    
    fig = figure();
    
    title(title_name);
   	set(fig,'Units','normalized') 
    set(fig, 'Position', [(0.05 + (degree_i-1) * 0.3) 0.25 0.3 0.5]);
    set(fig, 'DoubleBuffer', 'on');
    set(fig, 'Name', strrep(file_paths{degree_i},'_',' '));
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
    whitebg(fig,'black');
    set(fig, 'color', [0 0 0]);    
    
    hold on;
    
    if ~isempty(true_strands)
      
      x_components = [];
      y_components = [];
      z_components = [];
      
      for strand_i = 1:(size(true_strands,1)-1)
        
        true_strand = true_strands{strand_i,1};
        
        x_components = [x_components; true_strand(degree_i, 1)];
        y_components = [y_components; true_strand(degree_i, 2)];
        z_components = [z_components; true_strand(degree_i, 3)];        
      end
      
      plot3(x_components, y_components, z_components, 'o', 'MarkerEdgeColor','white', 'MarkerFaceColor','yellow', 'MarkerSize',8);

    end
    
    for strand_i = include_strands
   
      samples = strands{strand_i,1}; 

      num_samples = size(samples,1);
      
      if (isempty(include_samples))
        include_samples = [1:num_samples]';
        
      else
        include_samples = sort(include_samples(find(include_samples <= num_samples)))';
      end
      
      num_include_samples = length(include_samples);
      
      if invert_shades
        shades = ones(num_include_samples,1) - [end_shade:-(shade_range / (num_include_samples-1)):start_shade]';
        cdata = ones(num_include_samples, 3) - shades * ([1 1 1] - colours_of_evolved_strands(strand_i,:));
      else
        shades = [start_shade:(shade_range / (num_include_samples-1)):end_shade]';
        cdata = shades * (colours_of_evolved_strands(strand_i,:));
      end
      
      cdata = reshape(cdata, [num_include_samples 1 3]);
      
      h = surface(...
        'XData',[samples(include_samples,1) samples(include_samples,1)],...
        'YData',[samples(include_samples,2) samples(include_samples,2)],...
        'ZData',[samples(include_samples,3) samples(include_samples,3)],...
        'CData',[cdata cdata],...
        'FaceColor','none',...
        'EdgeColor','flat',...
        'Marker',marker, ...
        'LineStyle',style, ...        
        'LineWidth',2);
      view(3)     
      daspect ([ 1 1 1 ]);
    end
    
    hold off;
    

    daspect ([ 1 1 1 ]);
    
    cameratoolbar('Show');
    cameratoolbar('SetMode','orbit');
    
    
  end

  if display_colours
    display_colour_key(colours_of_evolved_strands);
  end  
    
  
end
