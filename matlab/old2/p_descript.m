function samples = p_descript(num_samples, coeff, varargin)
% samples = p_descript(num_samples, coeff, [plot_type])
% 
%	Args:			
%
%			num_samples - number of samples
%			coeff - the coefficients of the polynomial
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

range = ceil((num_samples - 1)/2);
num_samples = 2 * range + 1;

t = [-range:1:range]'./ range;




t_bar = ones(size(t));

if (degree >= 0)
    psi = ones(num_samples,1);
end


for (d = [1:(degree-1)])
	
    t_bar = t_bar .* t;
    
	psi = [psi, t_bar];

end

samples = psi * coeff;

%figure;
sfigure(1);
plot3(samples(:,1), samples(:,2), samples(:,3), plot_type);

cameratoolbar('Show');
cameratoolbar('SetMode','orbit');
grid;
xlabel('x');
ylabel('y');
zlabel('z');

max_lim = max(sum(coeff));

% xlim([-max_lim max_lim]);
% ylim([-max_lim max_lim]);
% zlim([-max_lim max_lim]);

% xlim([-100 100]);
% ylim([-100 100]);
% zlim([-100 100]);



coeff


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
