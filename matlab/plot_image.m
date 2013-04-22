function fig = plot_image(varargin)
%  
%  
% PURPOSE: Plots DW-MRI signals from a neighbourhood of voxels on a 3D grid
%  
%  
% ARGUMENTS: 
%  
%           image_filename     The filename of the generated image in .mif format.
%  
%  
% OPTIONS (name, description, type, default):
%  
%          -no_neg_values     No negative values will be plotted
%                              bool
%                              0
%  
%          -corner_offset     Offset of the bottom-left-back corner of the neighbourhood
%                              matrix_1x3
%                              minus half the neighbourhood dimension length
%  
%          -scale_signal      Scaling factor will be scaled by
%                              float
%                              automatically determined from signal
%  
%          -colour            Colour of the non-negative lobes of the signal
%                              matrix_1x3
%                              [1 1 0]
%  
%          -transparency      Transparency of the lobes
%                              float
%                              1
%  
%          -neg_colour        Colour of the negative lobes of the signal
%                              matrix_1x3
%                              [1 1 1]
%  
%          -grad_directions   Gradient encoding used.
%                              string
%                              '/data/home/tclose/Data/Gradient_directions/encoding_60.b'
%  
%          -plot_directions   The direction encoding used for the plotting.
%                              string
%                              '/data/home/tclose/Data/Gradient_directions/encoding_8000.b'
%  
%          -lmax              The lmax of the gradient encodings used.
%                              natural
%                              8
%  
%  



  description = 'Plots DW-MRI signals from a neighbourhood of voxels on a 3D grid';
  
  arguments = {'image_filename', 'The filename of the generated image in .mif format.'};
  
  if strfind('-help_display', arguments{1}) == 1

    auto_offset = 'auto_scale';
    auto_scale = 'auto_offset';
    auto_dim = 'auto_dim';
    auto_plot_directions ='depends on number of voxels';
    
    help_display = 1;
    
  else
    
    image_filename = varargin{1};

    img_struct = load_image(image_filename);

    if (~isfield(img_struct, 'data'))
      error(['Could not read image from file ' image_filename]);
    end 

    b0_mask = img_struct.DW_scheme(:,4) == 0.0;
    img = img_struct.data(:,:,:,~b0_mask);
    max_intensity = max(max(max(max(abs(img)))));
    max_b0 = max(max(max(max(abs(img_struct.data(:,:,:,b0_mask))))));
              
    if max_intensity == 0
      error('Neighbourhood contains no signal');
    end  

    if isfield(img_struct, 'offsets')
      auto_offset = img_struct.vox(1:3) .* str2num(img_struct.offsets);  %#ok<ST2NM>
    else
      auto_offset = [(-img_struct.vox(1) * img_struct.dim(1) * 0.5), (-img_struct.vox(2) * img_struct.dim(2) * 0.5), (-img_struct.vox(3) * img_struct.dim(3) * 0.5)];
    end
    auto_scale = min(img_struct.vox) / ( 2.25 * max_intensity);
    auto_dim = img_struct.dim(1:3);
    
    num_voxels = img_struct.dim(1) * img_struct.dim(2) * img_struct.dim(3);
    
    if num_voxels <= 125
      auto_plot_directions = '/home/tclose/Data/Tractography/diffusion/encoding/encoding_1000.b';
    else
      auto_plot_directions = '/home/tclose/Data/Tractography/diffusion/encoding/encoding_60.b';
    end

  end  
  
  options = {...
            'no_neg_values  ', 0,       'bool', 'No negative values will be plotted';...
            'corner_offset  ', auto_offset,...
                                      'matrix_1x3', 'Offset of the bottom-left-back corner of the neighbourhood (for matching with fibres that are defined in a particular region in space3';...
            'dim'           , auto_dim, 'matrix_1x3', 'The dimensions of a subset of the original image to plot';...
            'index_offset'   , [0 0 0], 'matrix_1x3', 'The index offset that the subset of image voxels are taken from';...
            'scale_signal   ', auto_scale,...
                                      'float', 'Scaling factor will be scaled by';...
            'colour         ',       [1 1 0],  'matrix_1x3', 'Colour of the non-negative lobes of the signal';...
            'transparency   ', 1.0,      'float', 'Transparency of the lobes';...
            'neg_colour     ',   [1 1 1],  'matrix_1x3', 'Colour of the negative lobes of the signal';...
            'grad_directions',        '/home/tclose/Data/Tractography/diffusion/encoding/encoding_60.b',...
                                      'string', 'Gradient encoding used.';...            
            'plot_directions',        auto_plot_directions,...
                                      'string', 'The direction encoding used for the plotting.';...
            'lmax           ',         8,        'natural', 'The lmax of the gradient encodings used.';...
            'no_vox_lines   ',         0,        'bool', 'Whether to use vox lines to plot.'
            'maxes_only'     ,         0,        'bool', 'Don''t attempt to plot the image and only print the max values'};
                                     
  
  parse_arguments      
  if (help_display) 
    return;
  end
  
  if any((dim + index_offset) > img_struct.dim(1:3))
    error(['Dims (' num2str(dim) ') plus offset (' num2str(index_offset) ') exceeds image dimensions (' num2str(img_struct.dim) ').'])
  end
  
  if maxes_only
      fprintf([image_filename ':\n    Max: ' num2str(max_intensity) '\n    b=0: '...
          num2str(max_b0) '\n'])
      return
  end
  
  plot_directions_mat = load(plot_directions);
  
  gradient_scheme = gen_scheme(img_struct.DW_scheme(:,1:3), lmax);
  plot_scheme = gen_scheme(plot_directions_mat(:,1:3), lmax);
     
  disp('');
  disp('-------------------------------');
  disp(' Image Properties: ');
  disp('-------------------------------');
  if isfield (img_struct, 'num_points') 
    disp(['  Number of points: ' num2str(img_struct.num_points)]);  
  elseif isfield (img_struct, 'num_segments') 
    disp(['  Number of segments: ' num2str(img_struct.num_segments)]);
  end
  if isfield (img_struct, 'num_strands')   
    disp(['  Number of strands per tract dimension: ' num2str(img_struct.num_strands) '']);  
  end
  if isfield (img_struct, 'state_location')
    disp(['  State location: ''' img_struct.state_location '''']);  
  end
  if isfield (img_struct, 'interpolation')   
    disp(['  Interpolation type: ''' img_struct.interpolation '''']);  
  end
  if isfield (img_struct, 'response_coeffs')  
    disp(['  Response coefficients: ' img_struct.response_coeffs]);      
  end
  if isfield (img_struct, 'software_version')   
    disp(['  Software version: ' img_struct.software_version ]);  
  end  
  if isfield (img_struct, 'datetime')
    disp(['  Generated: ' img_struct.datetime ]);  
  end  
  disp('-------------------------------');  
  
  
  
  fig = my_figure([image_filename ': (Max: ' num2str(max_intensity) ', b=0: '  num2str(max_b0) ')'], 1, 1, [1 1 1], 1);
  
%   fig = figure();
% 
%   set(fig,'Units','normalized') 
%   set(fig, 'Position', [0.3 0.25 0.4 0.5]);
%   set(fig, 'DoubleBuffer', 'on');
%   whitebg(fig, 'black');
%   set(fig, 'color', [0 0 0]);
%   daspect ([ 1 1 1 ]);
% 
%   cameratoolbar('Show');
%   cameratoolbar('SetMode','orbit');  
  
  set(gca, 'Xlim', [corner_offset(1) corner_offset(1) + img_struct.vox(1) * img_struct.dim(1)]);
  set(gca, 'Ylim', [corner_offset(2) corner_offset(2) + img_struct.vox(2) * img_struct.dim(2)]);
  set(gca, 'Zlim', [corner_offset(3) corner_offset(3) + img_struct.vox(3) * img_struct.dim(3)]);
  
  xlabel('x-axis');
  ylabel('y-axis');
  zlabel('z-axis');
  
  hold on;
  
  count = 0;
  
  
  for (z = (index_offset(3)+1):(dim(3) + index_offset(3)))

    for (y = (index_offset(2)+1):(dim(2) + index_offset(2)))

      for (x = (index_offset(1)+1):(dim(1) + index_offset(1)))
        
        voxel_signal = squeeze(img(x,y,z,:));
        
        voxel_signal_SH = amp2SH(voxel_signal, gradient_scheme);

        
        S = SH2amp (voxel_signal_SH, plot_scheme); % The interpolated voxel signal.

        
        S = S .* scale_signal;

        S2 = -S;

        S(find(S<0)) = 0;
        
        if no_neg_values
          S2 = zeros(size(S));
        else
          S2(find(S2<0)) = 0;
        end


        vox_offset = [(x-0.5) * img_struct.vox(1) (y-0.5) * img_struct.vox(2) (z-0.5) * img_struct.vox(3)] + corner_offset;

        
        vox_offset = repmat(vox_offset, [size(plot_scheme.vert,1) 1]);
        
        
        V = plot_scheme.vert .* (S*[ 1 1 1 ]) + vox_offset;
        V2 = plot_scheme.vert .* (S2*[ 1 1 1 ]) + vox_offset;


        h = patch('Vertices', V, 'Faces', plot_scheme.mesh); 
        set(h, 'LineStyle', 'None', ...
            'FaceLighting', 'Phong', ...
            'FaceColor', colour, ...
            'FaceAlpha', transparency);

        h2 = patch('Vertices', V2, 'Faces', plot_scheme.mesh); 
        set(h2, 'LineStyle', 'None', ...
                'FaceLighting', 'Phong', ...
                'FaceColor', neg_colour, ...
                'FaceAlpha', transparency);

        light('position', [1 1 1]);
        
        count = count + 1;
        
        fprintf('.');
        if mod(count,100) == 0
          fprintf('\n')
        end

      end
    end
  end
  
  hold off;  
  
   if ~no_vox_lines
    add_vox_lines_to_plot(img_struct.vox(1), max(dim));
  end

  lighting gouraud;

  fprintf('\n');  

  
end


