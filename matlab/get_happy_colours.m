function colours = get_happy_colours(colours, happy_colours, inv_happy_colours)
  if happy_colours
    if inv_happy_colours
      error('Can''t use ''-happy_colours'' and ''-inv_happy_colours'' simultaneously');
    end
    load([getuserdir(), '/Data/Tractography/misc/comb_happy_colours.mat']);
    colours = colours_of_bundles;
  elseif inv_happy_colours
    load([getuserdir(), '/Data/Tractography/misc/inv_comb_happy_colours.mat']);
    colours = colours_of_bundles;
  end
end

