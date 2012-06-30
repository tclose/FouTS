function batch_curve_prior_est(dirname)
% function batch_curve_prior_est(collection_name)
% 
% Estimates the prior distribution of Fourier coefficients for a given phantom 
% collection
% 
% 
%     dirname - directory name of the strand collection


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
				
                    prior_estimation(collection_dir);
					
					
				end
            end


        end    

    end
	
	
end