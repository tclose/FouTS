function save_samples_for_talk(output_path, camera_viewangle, paper_size, resolution, camera_target, camera_pos)

% 	set_transparent = 0;
 	
  if ~exist('resolution','var')
    resolution = 600;  
  end

  if ~exist('paper_size','var')
    paper_size = 15;
  end
  
  if ~exist('camera_pos','var')
    camera_pos = [-3.25   -3.25    3.25];
  end
  
  if ~exist('camera_viewangle','var')
    camera_viewangle = 7.5;
  end
  
  if ~exist('camera_target','var')
    camera_target = [0,0,0];
  end
  
  set(gca,'CameraPositionMode', 'manual');
  set(gca,'CameraTargetMode', 'manual');
  set(gca,'CameraViewAngleMode', 'manual');
    
  set(gca, 'CameraPosition', camera_pos);
  set(gca, 'CameraTarget', camera_target);
  set(gca, 'CameraViewAngle', camera_viewangle);

  
  %set(gcf, 'color', [1 1 1]);

  %note A4 is 21cm wide 29.7cm tall

  % xlim([-2 2]);
  % ylim([-2 2]);
  % zlim([-2 2]);

  %set(gcf,'Units','centimeters');
  set(gcf,'PaperUnits','centimeters');

  % set(gcf, 'Position', [10 10 10 10]);
  set(gcf, 'PaperSize', [paper_size paper_size]);
  set(gcf, 'PaperPosition', [0 0 paper_size paper_size]); 

  set(gca, 'visible', 'off');
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