function fourier_prior(collection_name,degree)
% function fourier_prior(dirname)
% 
% Estimates the prior distribution of Fourier coefficients for a given phantom 
% collection
% 
% 
%     collection_name - name of directory that holds the strand collections
%     degree - Degree of the Fourier descriptors


    [coeffs, residuals, num_control_points]= load_fourier_coeff(collection_name, degree);
    
    num_strands = size(coeffs,3);
    
    flat_coeffs = reshape(coeffs, degree, (3 * num_strands));
    	
%     for (degree_i = 1:degree)
%         
%        figure(degree_i);
%        
%        hist(flat_coeffs(degree_i,:)')
%        title(['Coeff ' num2str(degree_i)]);
%         
%     end
    
    figure(degree+1);
    hist(residuals);
    
    
    
end