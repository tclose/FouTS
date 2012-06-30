function samples = fr_descript(num_samples, coeff, varargin)
% samples = fr_descript(num_samples, coeff, [plot_type])
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


theta = 2 * pi * cos(k * pi) ;
phi = pi * cos(2 * k * pi);

samples = zeros(num_samples,3);

for (i=2:num_samples)
	
	samples(i,1) = (cos(theta(i)) + samples(i-1,1));
	samples(i,2) = (sin(theta(i)) + samples(i-1,2));
	samples(i,3) = (sin(phi(i)) + samples(i-1,3));
end



%figure;
plot3(samples(:,1), samples(:,2), samples(:,3), plot_type);

cameratoolbar('Show');
cameratoolbar('SetMode','orbit');

coeff
% xlim(plot_limits(1,:));
% ylim(plot_limits(2,:));
% zlim(plot_limits(3,:));