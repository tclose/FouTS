function centre_voxel = plot_centre_voxel(image_filename)

	load('/data/home/tclose/Code/Tractography/matlab_scripts/grad_scheme');
	load('/data/home/tclose/Code/Tractography/matlab_scripts/p_scheme');


	img_struct = read_image(image_filename);
	img = img_struct.data;

	centre_voxel = squeeze(img(10,10,10,:));
	sh = amp2SH(centre_voxel, grad_scheme);
	
	plot_SH(sh, p_scheme);

	daspect([1,1,1]);
	
end
