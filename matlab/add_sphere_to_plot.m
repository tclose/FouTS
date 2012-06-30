function add_sphere_to_plot(sphere_radius, transparency, props, num_divisions, centre)
  
  if ~exist('num_divisions', 'var')
    num_divisions = 100;
  end

  if ~exist('centre', 'var')
    centre = [0;0;0];
  end

  if ~exist('transparency','var')
    transparency = 0.25;
  end
  
  
  if sphere_radius == -1

    if exist('props', 'var')

      if isfield(props, 'prior_end_on_sphere_radius')
        sphere_radius = str2double(props.prior_end_on_sphere_radius);
      elseif isfield(props, 'prior_end_out_sphere_radius')
        sphere_radius = str2double(props.prior_end_out_sphere_radius);
      elseif isfield(props, 'prior_end_in_sphere_radius')
        sphere_radius = str2double(props.prior_end_in_sphere_radius);
      else
        %Do nothing.
        return;
      end
    end
  end


  if sphere_radius ~= 0
  
    num_faces = round(sphere_radius * num_divisions);

    [X, Y, Z] = sphere(num_faces);

    X = X * sphere_radius + centre(1);
    Y = Y * sphere_radius + centre(2);
    Z = Z * sphere_radius + centre(3);

    hold on;    

    h = surf(X, Y, Z);

    hold off;

    set(h,'facecolor', [0.5 0.5 0.5]);
    set(h,'FaceAlpha', transparency);
    set(h,'edgecolor', 'none');
  end

  
end