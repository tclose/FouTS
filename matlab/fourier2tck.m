function tck = fourier2tck(fourier, varargin)
% function tck = fourier2tck(fourier, [num_points=50])
% 
% converts fourier descriptors into a track of num_points length.
% 

if (nargin == 2)
	num_points = varargin{1};
else
	num_points = 100;
end

degree = size(fourier,1);

k = [0:1:(num_points-1)]'./(num_points-1);


if (degree < 1)
	error('Fourier coefficients have degree of less than 1');
end

psi = ones(num_points,1);

for (d = [1:(degree-1)])
	
	psi = [psi, sqrt(2) * cos(pi * k * d)];

end

tck = psi * fourier;