function params_evolution = load_parameter_evolution(dir_name)


	files = dir(dir_name)'; %'

	if (size(files,2) == 0) 
		error(['Could not load any strands from directory ' dirname ]);
  end

  sample_files = cell(0);
  
	for file = files

		if (~file.isdir)

			if (match_file_ext(file.name, 'frr' ))
        
        delimeters = [strfind(file.name, '_') strfind(file.name, '.frr')];

        if (length(delimeters) == 2 && strmatch(file.name(1:delimeters(1)), 'sample_'))
          sample_files{end+1} = file.name((delimeters(1)+1):(delimeters(2)-1)); %#ok<AGROW>
        end
        
      end
    end
  end
  
  sample_files = sample_files';

  num_samples = size(sample_files,1);
  
  max_decimal_places = 0;
  
  for file_i = 1:num_samples
    
    decimal_places = length(sample_files{file_i});
    
    if (decimal_places > max_decimal_places)
      
       max_decimal_places = decimal_places;
       
    end
  end
  
  for file_i = 1:num_sample_files
    
    decimal_places = length(sample_files{file_i});
    
    for i = 1:(max_decimal_places - decimal_places + 1)
    
      sample_files{file_i} = ['0' sample_files{file_i}];
      
    end
  end
  
  sample_files = sort(sample_files);

  fouriers = load_strands([dirname filesep 'sample_' sample_files(1) '.frr']); 
  
  num_strands = size(fouriers,1);
  degree = size(fouriers{1,1},1);  
  
  params_evolution = zeros(num_samples, num_strands, degree, 3);
  
  
  for sample_i = 1:num_samples

    fouriers = load_strands([dirname filesep 'sample_' sample_files(file_i) '.frr']); 
   
    for (strand_i = 1:num_strands)
          
      fourier = fouriers{strand_i,1};

      params_evolution(sample_i, strand_i, :, :) = reshape(fourier, [1 1 size(fourier,1) 3]);

    end
    
  end
  
  
end

