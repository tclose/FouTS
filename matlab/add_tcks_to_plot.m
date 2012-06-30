function add_tcks_to_plot(tcks, radii, colours, bundle_indices, num_tube_corners)

  if ~exist('num_tube_corners','var')
    num_tube_corners = 12;
  end
  
  hold on;
  
  for tck_i = 1:length(tcks)

    tck = tcks{tck_i};
    
    num_control_points = size(tck,1);

    segs = tck(2:num_control_points,:) - tck(1:num_control_points-1,:);

    norm_segs = zeros(size(segs));

    length_segs = sqrt(dot(segs,segs,2));

    zero_length_segs = find(length_segs == 0.0);
    
    if (size(zero_length_segs,1) >= 1) 
      length_segs(zero_length_segs) = 0.0001;
      segs(zero_length_segs) = 0.0001;
      tck(zero_length_segs) = tck(zero_length_segs)  + 0.0001; 
    end

    norm_segs(:,1) = segs(:,1) ./ length_segs;
    norm_segs(:,2) = segs(:,2) ./ length_segs;
    norm_segs(:,3) = segs(:,3) ./ length_segs;

    consec_segs_align = dot(norm_segs(1:num_control_points-2,:), norm_segs(2:num_control_points-1,:),2);

    count = 0;
    while (1)

      ref_vector = rand(1,3);					

      [t,n,b] = frame( tck(:,1),tck(:,2),tck(:,3), ref_vector);

      if (all(~isnan(n)))
        break;
      end

      if (count > 100) 
        error(['likely nan values in tck ' num2str(tck_i-1) '!']);
      end

      count = count + 1;
    end



    X = zeros(num_control_points, num_tube_corners);
    Y = zeros(num_control_points, num_tube_corners);
    Z = zeros(num_control_points, num_tube_corners);

    theta = 0:(2*pi/(num_tube_corners-1)):(2*pi);

    count = 0;

    for point_i = 1:num_control_points

      if (point_i == 1) 
        w = tck(1, :) + n(1,:);
        n_prime = n(1,:);
        b_prime = b(1,:);

      else

        mu = dot(t(point_i,:), tck(point_i,:) - w, 2) / dot(t(point_i,:), segs(point_i-1,:),2);

        w_proj = w + mu * segs(point_i-1, :);

        n_prime = w_proj - tck(point_i,:);

        n_prime = n_prime ./ norm(n_prime);

        b_prime = cross( t(point_i,:), n_prime);

        b_prime = b_prime ./ norm(b_prime);

        w = tck(point_i,:) + n_prime;

      end

      X(point_i,:) = tck(point_i, 1) + radii(tck_i) * ( n_prime(1,1) * cos(theta) + b_prime(1,1) * sin(theta));
      Y(point_i,:) = tck(point_i, 2) + radii(tck_i) * ( n_prime(1,2) * cos(theta) + b_prime(1,2) * sin(theta));
      Z(point_i,:) = tck(point_i, 3) + radii(tck_i) * ( n_prime(1,3) * cos(theta) + b_prime(1,3) * sin(theta));


    end

    h = surf(X,Y,Z);
    
    set(h,'facecolor', colours(bundle_indices(tck_i)+1,:));
    set(h,'edgecolor', 'none');

    if (mod(tck_i,25) == 0)
      fprintf('.');
    end

  end
  
  hold off;
  
end



function [t,n,b]=frame(x,y,z,vec)

% FRAME Calculate a Frenet-like frame for a polygonal space curve
% [t,n,b]=frame(x,y,z,v) returns the tangent unit vector, a normal
% and a binormal of the space curve x,y,z. The curve may be a row or
% column vector, the frame vectors are each row vectors. 
%
% This function calculates the normal by taking the cross product
% of the tangent with the vector v; if v is chosen so that it is
% always far from t the frame will not twist unduly.
% 
% If two points coincide, the previous tangent and normal will be used.
%
% Written by Anders Sandberg, asa@nada.kth.se, 2005


	N=size(x,1);
	if (N==1)
	  x=x'; %'
	  y=y'; %'
	  z=z'; %'
	  N=size(x,1);
	end

	t=zeros(N,3);
	b=zeros(N,3);
	n=zeros(N,3);

	p=[x y z];

	for i=2:(N-1)
	  t(i,:)=(p(i+1,:)-p(i-1,:));
	  tl=norm(t(i,:));
	  if (tl>0)
		t(i,:)=t(i,:)/tl;
	  else
		t(i,:)=t(i-1,:);
	  end
	end

	t(1,:)=p(2,:)-p(1,:);
	t(1,:)=t(1,:)/norm(t(1,:));

	t(N,:)=p(N,:)-p(N-1,:);
	t(N,:)=t(N,:)/norm(t(N,:));

	for i=2:(N-1)
	  n(i,:)=cross(t(i,:),vec);
	  nl=norm(n(i,:));
	  if (nl>0)
		n(i,:)=n(i,:)/nl;
	  else
		n(i,:)=n(i-1,:);
	  end
	end

	n(1,:)=cross(t(1,:),vec);
	n(1,:)=n(1,:)/norm(n(1,:));

	n(N,:)=cross(t(N,:),vec);
	n(N,:)=n(N,:)/norm(n(N,:));

	for i=1:N
	  b(i,:)=cross(t(i,:),n(i,:));
	  b(i,:)=b(i,:)/norm(b(i,:));
	end
end