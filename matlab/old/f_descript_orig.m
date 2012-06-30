function samples = f_descript_orig(num_samples, coeff)

if (size(coeff,2) ~= 3) 
	error(['The second dimension of the coefficient matrix must be 3, it was ' num2str(size(coeff,2))]);
end


k = [1:num_samples]'./num_samples;

psi = ones(num_samples,1);

for (degree = [1:(size(coeff,1)-1)])
	
	psi = [psi; sqrt(2) * cos (pi * degree *k)];
		
end


samples = psi * coeff;

plot3(samples(:,1), samples(:,2), samples(:,3), 'bx');

cameratoolbar('Show');
cameratoolbar('SetMode','orbit');
