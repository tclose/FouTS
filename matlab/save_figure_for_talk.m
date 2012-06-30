function save_figure_for_talk(output_path, resolution)

% 	set_transparent = 0;
 	
  if ~exist('resolution','var')
    resolution = 600;  
  end
 	 

  %set(gcf, 'color', [1 1 1]);

  %note A4 is 21cm wide 29.7cm tall

  % xlim([-2 2]);
  % ylim([-2 2]);
  % zlim([-2 2]);

  %set(gcf,'Units','centimeters');
  set(gcf,'PaperUnits','centimeters');

  % set(gcf, 'Position', [10 10 10 10]);
  set(gcf, 'PaperSize', [3 3]);
  set(gcf, 'PaperPosition', [0 0 3 3]); 

  %set(gca, 'visible', 'off');
  set(gca, 'Outerposition', [0 0 1 1]);
  set(gca, 'Position', [0 0 1 1]);

  background = [0 0 0];

  set(gcf, 'color', background);

  set(gcf,'InvertHardCopy','off');

  print(gcf, '-dpng',  ['-r' num2str(resolution)], [output_path]); 

  cdata = imread([output_path '.png']);
  %cdata = imcrop(cdata, [10 10 300 300]);

  % write it back out - setting transparency info
  imwrite(cdata, [output_path '.png'], 'png', 'BitDepth', 16, 'transparency', background) 
  
end