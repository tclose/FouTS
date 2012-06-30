function samples = f_descript(num_samples, coeff, varargin)
% samples = f_descript(num_samples, coeff, [plot_type])
% 
%	Args:			
%
%			samples - the samples, which are plotted
%			num_samples - number of samples
%			plot_type (optional - passed directly to plot3 (e.g. 'x' plots
%												crosses instead of a line)

if (nargin == 3) 
	plot_type = varargin{1};
else
	plot_type = '-';
end
			
if (size(coeff,2) ~= 3) 
	error('The second dimension of the coefficient matrix must be 3');
end

degree = size(coeff,1);

k = [0:1:(num_samples-1)]'./(num_samples-1);


if (degree >= 0)
    psi = ones(num_samples,1);
end

% start_point = coeff(1,:);
% end_point = coeff(1,:);

for (d = [1:(degree-1)])
	
	psi = [psi, sqrt(2) * cos(pi * k * d)];

% 	start_point = start_point + coeff(d,:);
% 	end_point = end_point + coeff(d,:) * cos(pi * d);
end

samples = psi * coeff;

%sfigure;
plot3(samples(:,1), samples(:,2), samples(:,3), plot_type);

%  = (1,:)' - 1.5 * scalars(2);
% up_lim = c(1,:)' + 1.5 * scalars(2);
% 
% xlim([low_lim(1) up_lim(1)]);
% ylim([low_lim(2) up_lim(2)]);
% zlim([low_lim(3) up_lim(3)]);    

cameratoolbar('Show');
cameratoolbar('SetMode','orbit');
grid;
xlabel('x');
ylabel('y');
zlabel('z');
coeff
% xlim(plot_limits(1,:));
% ylim(plot_limits(2,:));
% zlim(plot_limits(3,:));


function h = sfigure(h)
% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%
% See also figure

if nargin>=1 
	if ishandle(h)
		set(0, 'CurrentFigure', h);
	else
		h = figure(h);
	end
else
	h = figure;
end