function params_evolution = load_param_evolution(dir_name, strand_i, degree_i, dim_i)


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
  
  
  sample_files = sort(sample_files)';
  
  num_samples = size(sample_files,1);
  
  path = [ dir_name filesep 'sample_' sample_files{1} '.frr' ];
  
  path
  
  fouriers = load_strands(path); 
  
  num_strands = size(fouriers,1);
  degree = size(fouriers{1,1},1);  
  
  params_evolution = zeros(num_samples);
  
    
  for sample_i = 1:num_samples
    
    path = [dir_name, filesep, 'sample_', sample_files{sample_i}, '.frr'];
    
    fouriers = load_strands(path); 

    fourier = fouriers{strand_i};
    
    params_evolution(sample_i) = fourier(degree_i, dim_i);


  end

  
  plot(params_evolution);
  
  
end

