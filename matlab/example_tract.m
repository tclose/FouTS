function example_tract


  function prep_figure_for_save(cam_orbit, limits)
    
    xlim(limits(1,:));
    ylim(limits(2,:));
    zlim(limits(3,:));
    
    daspect ([ 1 1 1 ]);
    camorbit(cam_orbit(1),cam_orbit(2));

    cameratoolbar('Show');
    cameratoolbar('SetMode','orbit');

    light           
    lighting gouraud;

  end

  close all;

  
%   Set defaults
  
  red_colour = [ 0.7, 0.2, 0.2];
  grey_colour = [0.7,0.7,0.7];
  green_colour = [0.2,0.7,0.2];
  blue_colour = [0.2,0.2,0.7];

  num_width_sections = 6;

  strand_radii = 0.005;

  cam_orbit = [70, -65];
  
  def_dlimit = .275;
  
  data_limits = [-def_dlimit, def_dlimit; -def_dlimit, def_dlimit; -def_dlimit, def_dlimit];

  
  
%   Load data
  
  [tract, prop_keys, prop_value]    = load_tracts('/home/tclose/Data/Tractography/poster/example/deform.tct');
  base_width = get_properties(prop_keys, prop_value, 'base_width', 0.01);  

  backbone = load_strands('/home/tclose/Data/Tractography/poster/example/backbone.str');
  back_to_strand1 = load_strand_sets('poster/example/back_to_strand1.sst');
  back_to_strand2 = load_strand_sets('poster/example/back_to_strand2.sst');
  strand1_about_back = load_strand_sets('poster/example/str1_about_back.sst');
  strand2_about_back = load_strand_sets('poster/example/str2_about_back.sst');
  

% Backbone
  
  back_fig = my_figure('Backbone');

  add_strands_to_plot(backbone, strand_radii, red_colour , [0]);

  prep_figure_for_save(cam_orbit, data_limits);

  save_figure_for_poster('/home/tclose/Documents/Tractography/ismrm/images/backbone');

  
 
 % Backbone to Strand1
   
   back_to_str1_fig = my_figure('Backbone to Strand1');
 
   add_strands_to_plot(backbone, strand_radii, red_colour , [0]);
   
   for set_i = 2:size(back_to_strand1,1)
     add_strands_to_plot(back_to_strand1{set_i}, strand_radii, grey_colour , [0]);
   end
 
   prep_figure_for_save(cam_orbit, data_limits);
 
   save_figure_for_poster('/home/tclose/Documents/Tractography/ismrm/images/back_to_str1');
 
   
   
 % Backbone to Strand2
 
   back_to_str2_fig = my_figure('Backbone to Strand2');
 
   
   add_strands_to_plot(backbone, strand_radii, red_colour , [0]);
 
   
   for set_i = 2:size(back_to_strand2,1)
     add_strands_to_plot(back_to_strand2{set_i}, strand_radii, grey_colour , [0]);
   end
   
   prep_figure_for_save(cam_orbit, data_limits);
 
 
   save_figure_for_poster('/home/tclose/Documents/Tractography/ismrm/images/back_to_str2');
 

    

% Strand1 about Backbone 
  
  num_strands = size(strand1_about_back,1);
  
  half = floor(num_strands / 2);
  
  str1_about_back_fig = my_figure('Strand1 about Backbone');
  
  add_strands_to_plot(backbone, strand_radii, red_colour , [0]);
  
  for set_i = [1:half, (half+2):num_strands]
    add_strands_to_plot(strand2_about_back{set_i}, strand_radii, grey_colour , [0]);
  end

  prep_figure_for_save(cam_orbit, data_limits);

  save_figure_for_poster('/home/tclose/Documents/Tractography/ismrm/images/str1_about_back');

  
  
% Strand2 about Backbone

  num_strands = size(strand2_about_back,1);
  
  half = floor(num_strands / 2);

  str2_about_back_fig = my_figure('Strand2 about Backbone'); 

  add_strands_to_plot(backbone, strand_radii, red_colour , [0]);
  
  for set_i = [1:half, (half+2):num_strands]
    add_strands_to_plot(strand1_about_back{set_i}, strand_radii, grey_colour , [0]);
  end
  
  prep_figure_for_save(cam_orbit, data_limits);


  save_figure_for_poster('/home/tclose/Documents/Tractography/ismrm/images/str2_about_back');

  
  
  
% Tract Strands



  [strands, bundle_indices] = tracts2strands(tract, base_width, num_width_sections);

  num_strands = length(strands);
  
  bundle_colours = repmat(grey_colour, num_strands,1);
  
  central_strand_index = (num_strands-1)/2+1;
  
  
  for strand_i = 31:39
    bundle_colours(strand_i,:) = blue_colour;
  end

  for strand_i = [3,9,17,26,44,53,61,67]
    bundle_colours(strand_i,:) = green_colour;
  end
  
  bundle_colours(central_strand_index,:) = red_colour;
  

%   set_bundle_colours(num_strands);
%  
%   global colours_of_bundles;
%   bundle_colours = colours_of_bundles;

  
  tract_fig = my_figure('Tract Strands');


  add_strands_to_plot(strands, repmat(strand_radii, num_strands,1),  bundle_colours, [0:1:num_strands-1]);  


  prep_figure_for_save(cam_orbit, data_limits);


  save_figure_for_poster('/home/tclose/Documents/Tractography/ismrm/images/tract_strands');

  
  
% Whole Tract 

  tract_fig = my_figure('Whole Tract');

  add_tracts_to_plot(tract, grey_colour, 1.0, base_width,0,100,64);

  prep_figure_for_save(cam_orbit, data_limits);

  save_figure_for_poster('/home/tclose/Documents/Tractography/ismrm/images/whole_tract');
  
end