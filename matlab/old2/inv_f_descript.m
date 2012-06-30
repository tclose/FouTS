function coeff = inv_f_descript(samples, degree, varargin)
% samples = inv_f_descript(samples, degree, [plot_type])
% 
%	Args:			
%
%			coeff - the fourier coefficients
%			samples - the samples of the line to be transformed
%			degree - the degree of the basis set
%			plot_type (optional - passed directly to plot3 (e.g. 'x' plots
%												crosses instead of a line)

if (size(samples,2) ~= 3) 
	error('The second dimension of the samples matrix must be 3');
end

num_samples = size(samples, 1);

if (nargin >= 3) 
	plot_type = varargin{1};
else
	plot_type = '-';
end

if (nargin == 4)
	num_output_samples = varargin{2};
else
	if (num_samples > 100)
		num_output_samples = num_samples;
	else
		num_output_samples = 100;
	end
end



k = [0:1:(num_samples-1)]'./(num_samples-1);

psi = ones(num_samples,1);


for (d = [1:(degree-1)])
    
    psi = [psi, sqrt(2) * cos(pi * k * d)];
        
end
    
coeff = inv(psi' * psi) * psi' * samples;

if (num_samples ~= num_output_samples)
	psi = [];
	
	k = [0:1:(num_output_samples-1)]'./(num_output_samples-1);
		
	psi = ones(num_output_samples,1);
	
	for (d = [1:(degree-1)])

		psi = [psi, sqrt(2) * cos(pi * k * d)];

	end
	
end
	
new_samples = psi * coeff;

figure();
hold on;

plot3(samples(:,1),samples(:,2),samples(:,3), [plot_type 'b']); 

plot3(new_samples(:,1),new_samples(:,2),new_samples(:,3), [plot_type 'g']); 

hold off;


cameratoolbar('Show');
cameratoolbar('SetMode','orbit');