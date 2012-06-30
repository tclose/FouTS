function q = get_response_coeffs(SH)
%Takes a vector of spherical harmonic coefficients (SH) of a given signal model
%of diffusion and calculates the corresponding coefficients for the polynomial
%expansion used in the calculation of the instantaneous signal intensity.

    R = get_m0_coeffs(SH)

    max_order = 2 * (length(R)-1);

    l = [0:2:max_order]; 

    if (max_order > 8) 
        error(['Sorry only polynomials up to Legendre polynomial order 8 have been calculated (order ' num2str(max_order) ' supplied)']);
    end    

    P = zeros(5);

    P(1,:) =         [ 1     0    0      0    0];	    %0th order Associated Legendre polynomial, m=0
    P(2,:) = 1/2   * [-1     3    0      0    0];     %2nd order Associated Legendre polynomial, m=0
    P(3,:) = 1/8   * [ 3   -30   35      0    0];     %4th order Associated Legendre polynomial, m=0
    P(4,:) = 1/16  * [-5   105 -315    231    0];     %6th order Associated Legendre polynomial, m=0
    P(5,:) = 1/128 * [35 -1260 6930 -12012 6435];     %8th order Associated Legendre polynomial, m=0


    %Weight the Legendre polynomials by appropriate scaling function for
    %spherical harmonics. 
	
    
    normalise = diag(sqrt(  (2 * l + 1) ./ (4 * pi)  ));
    
    P = normalise * P;				

    P = diag(R) * P

    q = sum(P,1)';



end


function m0 = get_m0_coeffs(SH)

    reshape(SH, length(SH),1);


    total_index = 0;
    i = 0;
    indices = [];

    while (total_index < length(SH))   


        block_size = 2 * i;

        new_index = total_index + block_size + 1;


        indices = [indices; new_index];

        total_index = new_index + block_size;
        i = i + 1;
    end

    m0 = SH(indices);
end