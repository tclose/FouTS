function add_isotropic_region_to_plot(isotropic_regions)

	for (isotropic_region_i = 1:size(isotropic_regions,1))
		
		centre = isotropic_regions(isotropic_region_i, 1:3)';
		radius = isotropic_regions(isotropic_region_i, 4);
		
		num_faces = round(radius * 100);
		
		[X, Y, Z] = sphere(num_faces);
	
		X = X * radius + centre(1);
		Y = Y * radius + centre(2);
		Z = Z * radius + centre(3);
		
		h = surf(X, Y, Z);
		
		set(h,'facecolor', [0.5 0.5 0.5]);
		set(h, 'FaceAlpha', 0.5);
	end	



end