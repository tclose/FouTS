function coeffs = rand_f_descript(num_samples, iterations, scalars, varargin)
%rand_f_descript(num_samples, iterations, scalars, [normalise], [plot_type])
%
%
%   rand_f_descript generates multiple plots of fourier descriptors in
%   succession to give a feel for the types of curves that can be produced.
%
%   The function will randomly generate a number (specified by the input
%   'iterations') of different curves and plot them in sucession (press any
%   key to iterate to next plot). 
%
%   Coefficients are initially generated between 0 and 1 and are scaled via 
%   the 'scalars' input vector.  Each dimension is scaled identically. If 
%   the optional argument normalise is set to 1 then the second degree 
%   coefficients (the linear length) are normalised before they are scaled.
%
%             num_samples - number of samples to plot along each curve
%             iterations - number of random generations 
%             scalars - scalars of the relative coefficient degrees
%             normalise - if set to one, the coefficients are normalised 
%                           across each degree before they are scaled.
%             plot_type - passed directly to the plot function. Use 'x' to
%                           see placement of samples.


if (nargin >=4)
    normalise = varargin{1};
else
    normalise = 1;
end

if (nargin == 5)
    plot_type = varargin{2};
else
    plot_type = '-';
end

if (size(scalars,1) < size(scalars,2))
	scalars = scalars';
end

degree = size(scalars,1);

coeffs = [];

for (i =1:iterations)
    c = 2*(rand(degree, 3) - 0.5);
   
    length = sqrt(sum(c(2,:) * c(2,:)'));
    
    if (normalise)    
        for (d = 1:degree)
           row_length = sqrt(sum(c(d,:) * c(d,:)'));        
           c(d,:) = c(d,:) ./row_length;
        end
    else
       c(2,:) = c(2,:) ./ length;
    end
    
    c = [c(:,1) .* scalars, c(:,2) .* scalars, c(:,3) .* scalars];
    i
	
	sfigure(1);
    f_descript(num_samples, c, plot_type);
%     title('Strand path');
% 
% 	if (degree == 4) 
% 		figure(2);
% 		plotv(c);
% 		title('Coefficient vectors');
% 	end
	
	
    low_lim = c(1,:)' - 1.5 * scalars(2);
    up_lim = c(1,:)' + 1.5 * scalars(2);
    
    xlim([low_lim(1) up_lim(1)]);
    ylim([low_lim(2) up_lim(2)]);
    zlim([low_lim(3) up_lim(3)]);    

	coeffs = cat(3, coeffs,c);
	
    pause;

end



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
