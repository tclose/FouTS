function plot_amp(S, scheme, colour, transparency, ax)

% function plot_amp(S, scheme, transparency)
%
% Plot amplitudes 'S' for the set of directions given 
% in 'scheme', with optional colour and transparency
% specified.

 
% %set(gcf, 'DefaultAxesCameraViewAngleMode', 'manual', ...
%          'DefaultAxesCameraViewAngle', 20, ...
%          'DefaultAxesCameraTarget', [ 0 0 0 ], ...
%          'DefaultAxesCameraPosition', 6*range*[0 0 1], ...
%          'DefaultAxesVisible', 'on', ...
%          'DefaultAxesBox', 'on', ...
%          'DefaultAxesXLim', [ -range range ], ...
%          'DefaultAxesYLim', [ -range range ], ...
%          'DefaultAxesZLim', [ -range range ]);


if ~exist('colour')
  colour = [ 1 1 0 ];
end

if ~exist('transparency') 
  transparency = 1;
end

if ~exist('ax')
	main_fig = figure();
	
	set(main_fig,'Units','normalized') 

	set(main_fig, 'Position', [0.3 0.25 0.4 0.5]);
	set(main_fig, 'DoubleBuffer', 'on');
	set(main_fig, 'Name', 'Amplitude plot');

	cameratoolbar('Show');
	cameratoolbar('SetMode','orbit');
	clf
	set(gcf, 'Color', [1 1 1])
else
	subplot(ax);
	
end

range = 1.5*max(S);
negcolour = [ 1 1 1 ];

S2 = -S;
S(find(S<0)) = 0;
S2(find(S2<0)) = 0;

V = scheme.vert .* (S*[ 1 1 1 ]);
V2 = scheme.vert .* (S2*[ 1 1 1 ]);


h = patch('Vertices', V, 'Faces', scheme.mesh); 
set(h, 'LineStyle', 'None', ...
    'FaceLighting', 'Phong', ...
    'FaceColor', colour, ...
    'FaceAlpha', transparency);

h2 = patch('Vertices', V2, 'Faces', scheme.mesh); 
set(h2, 'LineStyle', 'None', ...
        'FaceLighting', 'Phong', ...
        'FaceColor', negcolour, ...
        'FaceAlpha', transparency);

light('position', [1 1 1]);

ylabel('y');
xlabel('x');
zlabel('z');

if (~exist('ax'))
	drawnow
end
