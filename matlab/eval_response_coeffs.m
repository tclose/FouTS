function S = eval_response_coeffs(response_coeff, scheme, eig_v)

%Normalise and ensure that eig_v is 3x1
eig_v = reshape((eig_v ./ norm(eig_v)),3,1);

%Convert gradient scheme to cartesian coordinates
grad_C = s2c([ scheme.el scheme.az 1+0*scheme.az ])';

%Tile eig_v to match gradient scheme size
eig_v = repmat(eig_v , 1, size(grad_C,2));

%Get the cosine(angle) between the gradients and the principle eigenvector
cos_angle = dot(grad_C, eig_v)';

%Evaluate the signal functions from the response coefficients and the
%cosine angles
S = response_coeff(1) + response_coeff(2) * cos_angle.^2 + ...
    response_coeff(3) * cos_angle.^4 + response_coeff(4) * cos_angle.^6+...
    response_coeff(5) * cos_angle.^8;

