function batch_fourier_est(dirname,degree)
% function batch_fourier_est(dirname)
% 
% Estimates the prior distribution of Fourier coefficients for a given phantom 
% collection
% 
% 
%     dirname - name of directory that holds the strand collections
%     degree - Degree of the Fourier descriptors


    collection_names = dir(dirname);

    stats = [];
    
    if (size(collection_names,1) < 1)
        error(['Did not read any collection names from ' dirname]);
    else
    
        for (coll_i =1:size(collection_names,1))
            
            collection_name = collection_names(coll_i).name;
            
            if (collection_name(1,1) ~= '.')
				
                collection_dir = [dirname '/' collection_name];
				
				if (isdir(collection_dir)) 
				
                    load_fourier_coeff(collection_dir,degree);
					
					
				end
            end


        end    

    end
	
	
end