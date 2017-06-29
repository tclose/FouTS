function add_tracts_to_plot (tracts, colours, intensities, base_widths, tube_corners, num_length_sections, fixed_transparency, bundle_indices, colour_indices)


  if ~exist('tube_corners', 'var')
    tube_corners = 12;
  end
  
  if tube_corners < 1
    error (['''tube_corners'' option must be greater than 1 (' num2str(tube_corners) ').']);
  end

  if ~exist('num_length_sections', 'var')
    num_length_sections = 100;
  end  
    
  if ~exist('fixed_transparency', 'var')
    fixed_transparency = 1.0;
  end

  if ~exist('bundle_indices', 'var')
      bundle_indices = [];
  end
  
  if isempty(bundle_indices)
      bundle_indices = 0:1:(size(tracts,1)-1);
  end

  if ~exist('colour_indices', 'var')
      colour_indices = [];
  end

  if isempty(colour_indices)
      colour_indices = 1:max(bundle_indices+1);
  end  
  
  if nargin > 9
    error(['Incorrect number of arguments (' num2str(nargin) '). Cannot have more than 8']);
  end

  num_tracts = size(tracts,1);    

  hold on;

  for tract_i=1:num_tracts

    if base_widths(tract_i) == 0
      disp(['Warning!! Base width of tract ' num2str(tract_i) ' is zero.']);
    end
    
    X = zeros(num_length_sections, tube_corners);
    Y = zeros(num_length_sections, tube_corners);
    Z = zeros(num_length_sections, tube_corners);

    theta = 0:(2*pi/(tube_corners-1)):(2*pi);

    backbone  = fourier2tck(tracts{tract_i,1}, num_length_sections);
    axis1     = fourier2tck(tracts{tract_i,2}, num_length_sections);
    axis2     = fourier2tck(tracts{tract_i,3}, num_length_sections);      

    count = 0;

    for point_i = 1:num_length_sections

      X(point_i,:) = backbone(point_i, 1) + base_widths(tract_i) * ( axis1(point_i,1) * cos(theta) + axis2(point_i,1) * sin(theta));
      Y(point_i,:) = backbone(point_i, 2) + base_widths(tract_i) * ( axis1(point_i,2) * cos(theta) + axis2(point_i,2) * sin(theta));
      Z(point_i,:) = backbone(point_i, 3) + base_widths(tract_i) * ( axis1(point_i,3) * cos(theta) + axis2(point_i,3) * sin(theta));

    end

    h = surf(X,Y,Z);

    set(h,'facecolor', colours(colour_indices(bundle_indices(tract_i)+1)+1,:));
    set(h,'edgecolor', 'none');
    
    if fixed_transparency ~= 0
      set(h, 'FaceAlpha', fixed_transparency);
    else
      set(h,'FaceAlpha', intensities(tract_i));
    end
    
  end

  hold off;

end