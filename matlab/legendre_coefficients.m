function get_response_coeffs(R)
%Takes a vector of spherical harmonic coefficients (R) from a given response
%function and calculates the corresponding coefficients for the polynomial
%expansion used in the calculation of the instantaneous signal intensity.

R = [1 1 1 1 1];

if (length(R) > 5) 
	error('Sorry only polynomials up to order 4 have been calculated.');
end

p{1} = [0 0 0 0 1];								%0th order Legendre polynomial
p{2} = 1/2 * [0 0 0 3 -1];						%2nd order Legendre polynomial
p{3} = 1/8 * [0 0 35 -30 3];					%4th order Legendre polynomial
p{4} = 1/16 * [0 231 -315 105 -5];				%6th order Legendre polynomial
p{5} = 1/128 * [6435 -12012 6930 -1260 35];		%8th order Legendre polynomial


for (i = 1:size(p,1))
	P = [P;p{1}];
end
%[p0;p2;p4;p6;p8];

P = P * sqrt((4 * size(P,1)+1)/2);

P = P * diag(R);

q = sum(P)'



