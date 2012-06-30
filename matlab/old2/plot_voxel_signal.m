function voxel_signal = plot_voxel_signal(image_filename, x, y, z)

	load('/data/home/tclose/Code/Tractography/matlab_scripts/grad_scheme');
	load('/data/home/tclose/Code/Tractography/matlab_scripts/p_scheme');


	img_struct = read_image(image_filename);
	img = img_struct.data;

	voxel_signal = squeeze(img(x,y,z,:));
	sh = amp2SH(voxel_signal, grad_scheme);
	
	plot_SH(sh, p_scheme);

	daspect([1,1,1]);
	
end
