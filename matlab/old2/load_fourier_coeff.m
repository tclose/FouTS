function [coeffs, residual_vars, num_strand_cpoints, length] = load_fourier_coeff(collection_name, degree)
% function [coeffs, residual_vars, num_strand_control_points, straight_line_length] = load_fourier_coeff(collection_name, degree)
% 
% Estimates the prior distribution of Fourier coefficients for a given phantom 
% collection
% 
% 
%       coeffs                      - A degree X 3 X num_strands matrix containing the loaded coefficients
%       residual_vars               - Residual variance after least squares fit 
%       num_strand_cpoints   - Number of control points per strand
%       straight_line_length        - Straight-line length between strand
%                                       endpoints
% 

%       collection_name - directory name of the strand collection
%       degree          - Degree of the Fourier descriptors


    collection = load_strands(collection_name);

    num_strands = size(collection,1);

    coeffs = [];
    residual_vars = [];
    num_strand_cpoints = [];
    length = [];
    
    for  strand_i = 1:num_strands
        
        strand = collection{strand_i,1};
        num_cpoints = size(strand,1);
        
        if (num_cpoints > degree) 

            k = [0:1:(num_cpoints-1)]'./(num_cpoints-1);

            psi = ones(num_cpoints,1);

            for d = 1:(degree-1)

                psi = [psi, sqrt(2) * cos(pi * k * d)];

            end
            
            C = inv(psi' * psi) * psi' * strand;
            
            residual = strand - psi * C;
            
            residual_var = sum(sum(residual.^2))/num_cpoints; 
                        
            coeffs = cat(3,coeffs,C);
            
            residual_vars = [residual_vars; residual_var]; 
            num_strand_cpoints  = [num_strand_cpoints; num_cpoints];            
            length = [length; norm(strand(num_cpoints,:) - strand(1,:))]; 
            
        end
        
    end
end    


function strand_collection = load_strands(dirname)

%  function strand_collection = load_strands(dirname)
% 
%   load_strand_collection.m
%   Numerical Fibre Generator
% 
%   Created by Tom Close on 19/02/08.
%   Copyright 2008 Tom Close.
%   Distributed under the GNU General Public Licence.
% 
% 
% 
%   This file is part of 'Numerical Fibre Generator'.
% 
%   'Numerical Fibre Generator' is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   'Numerical Fibre Generator' is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with 'Numerical Fibre Generator'.  If not, see <http://www.gnu.org/licenses/>
% 
% 
% 


	num_strands = 0;

	files = dir(dirname)'; %'

	if (size(files) == [1,0]) 
		error(['Could not load any strands from directory ' dirname ]);
	end

	for file = files

		if (~file.isdir)

			delimeters = [strfind(file.name, '_') strfind(file.name, '-') strfind(file.name, '.txt' )];

			if (length(delimeters) == 4 && strmatch('strand', file.name))

				num_strands = num_strands + 1;


				strand_collection{num_strands, 1} = load([dirname filesep file.name]);
				strand_collection{num_strands, 2} = str2num(file.name(delimeters(1) + 1 :delimeters(2) -1 ));
				strand_collection{num_strands, 3} = str2num(file.name(delimeters(3) + 2 :delimeters(4) -1 ));
				strand_collection{num_strands, 4} = str2num(file.name(delimeters(2) + 1 :delimeters(3) -1 ));


			end
		end

	end
end



