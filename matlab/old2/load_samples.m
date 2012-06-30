function [samples, ext_prop_keys, ext_prop_values, initial_set, true_set, sphere_radius] = load_samples(filename, varargin)
% samples = f_descript(num_samples, coeff, [plot_type])
% 
%	Args:			
%
%			samples - the samples, which are plotted
%			

	if (nargin > 2)
		error(['Incorrect number of arguments ' nargin ' expecting no more than 2']);
  end
	
  if (nargin == 2)
		sample_indices = varargin{1};
	else
		sample_indices = [];
  end
    
	
	dots_i = findstr(filename,'.');
	ext = filename(dots_i(end):end);
	
	if ~strcmp(ext, '.tsp')
    error(['File has extension ' ext ', required to be ''tsp'' for tract set samples file.']);
  end
		
  
  if any(sample_indices < 1)
    error(['Sample indices must be positive integers (min index: ' num2str(min(sample_indices)) ')']);
  end
  
	
	file = fopen(filename,'r');	
	
	[fname, mode, machine_format] = fopen(file);

	
	if (file == -1) 
		error([ 'Could not open file ' filename '!' ]);
	end

	line = fgetl(file);
	
	if (~strcmp(line(1:13),'mrtrix tracks'))
		error(['File, ' filename ' was not a valid mrtrix tracks file']);
	end
	
	offset = 0;
  num_points_per_set = 0;
  num_tracts = 0;
  degree = 0;
  num_samples = 0;
	true_set_provided = false;	
  sphere_radius = 0.39;
  
	for (line_i = 1:1000)	
		
    line =  fgetl(file);
    
    if strcmp(line(1:3), 'END')
      break;
    end
    
    if length(line) >= 4 && strcmp(line(1:4),'file')
      offset = str2double(line(9:end));
      
    elseif length(line) >= 18 && strcmp(line(1:18), 'num_points_per_set')
      num_points_per_set = str2double(line(21:end)); 
      
    elseif length(line) >= 8 && strcmp(line(1:8), 'set_size')
      num_points_per_set = str2double(line(11:end));           
    
    elseif length(line) >= 11 && strcmp(line(1:11), 'num_tracts')
      num_tracts = str2double(line(14:end));

    elseif length(line) >= 5 && strcmp(line(1:5), 'count')
      num_samples = str2double(line(8:end));      
      
    elseif length(line) >= 17 && strcmp(line(1:17), 'true_set_provided')
      true_set_provided = true;
      
 		elseif length(line) >= 22 && strcmp(line(1:22),'proposal_sphere_radius')
 			sphere_radius = str2num(line(25:end));    

 		elseif length(line) >= 8 && strcmp(line(1:8),'datatype')
 			datatype = line(11:end);    
      
      if (strcmp(datatype, 'Float32BE'))
        machine_format = 'ieee-be';
      elseif (strcmp(datatype, 'Float32LE'))
        machine_format = 'ieee-le';
      else
        error(['Unrecognised data format ''' datatype ''', should be Float32BE or Float32LE.']);
      end
      
    end
    
  end	
	
	if (~offset || ~num_points_per_set)
		error(['All required properties, ''file'' and ''num_points_per_set'', were not found after 1000 lines']);
  end
  
  

  if any(sample_indices > num_samples)
    error(['Sample index (' num2str(max(sample_indices)) ') out of range (' num2str(num_samples) ').']);
  end
  
	fseek(file, offset, 'bof');

  samples = cell(0,1);
  
  %Read true tract set if provided.
  if (true_set_provided)  
    true_set = read_tract_set(file, num_points_per_set, machine_format);
  else  
    true_set = [];
  end
  
  %Read initial tract set.
  initial_set = read_tract_set(file, num_points_per_set, machine_format);
    
  
  if isempty(initial_set)
    error(['No initial set was found in samples file']);
  end
    
  %Load extended properties
  [ext_prop_keys, ext_prop_values] = load_extend_elem_props([filename 'x']);
  
  if isempty(sample_indices)
    
    end_of_file = false;
    
    while ~end_of_file      
    
      samples{end+1,1} = read_tract_set(file, num_points_per_set, machine_format);
    
      end_of_file = isempty(samples{end});
      
    end
    
  else

    set_byte_size = num_points_per_set * 3 * 4; % times 3 for a point containing 3 numbers, and times 4 for the size of a float.

    prev_index = 0;

    for set_i = 1:length(sample_indices)

      sample_index = sample_indices(set_i);

      if (sample_index - prev_index ~= 1)
        if fseek(file, offset + sample_index * set_byte_size, 'bof') < 0
          error(['Sample index (' num2str(sample_index) ') is out of range of the file']);
        end
      end

      samples{set_i+2,1} = read_tract_set(file, num_tracts, num_points_per_set, machine_format);

      if end_of_file
        error(['End of file reached while attempting to load tract sample ' num2str(sample_indices(set_i)) '.']);
      end
      
      prev_index = sample_index;

    end
    
    ext_prop_values = ext_prop_values(sample_indices,:);
    
  end
  
end


function tract_set = read_tract_set(file, num_points_per_set, machine_format)

  block_size = num_points_per_set * 3; 

  [block, read_count] = fread(file, block_size, 'float32=>double', machine_format);

  if isinf(block(1)) || feof(file)
    tract_set = [];
    return;
  end  

  if read_count ~= block_size
    error(['''num_points_per_set'' (' num2str(num_points_per_set) ') property does not divide evenly into the size of file']);
  end

  block = reshape(block, 3, num_points_per_set)';
  
  tract_seperators = find(isinf(block(:,1)));
  
  num_tracts = size(tract_seperators,1);

  if tract_seperators(end) ~= num_points_per_set
    error(['''num_points_per_set'' (' num2str(num_points_per_set) ') does not align with tract seperators ([-inf, -inf, -inf]).']);
  end
  
  
  tract_set = cell(num_tracts,3);
  
  tract_start_index = 1;
  
  for tract_i = 1:num_tracts

    tract_end_index = tract_seperators(tract_i);
    
    tract = block(tract_start_index:(tract_end_index-1),:);

    axes_seperators = find(isnan(tract(:,1)));
    
    if length(axes_seperators) ~= 3
      error (['Expected 3 axes in tract ' num2str(tract_i) ', only found ' length(axes_seperators) '.']);
    end
    
    axes_start_index = 1;
    
    for axes_i = 1:3

      axes_end_index = axes_seperators(axes_i);

      tract_set{tract_i, axes_i} = tract(axes_start_index:(axes_end_index-1),:);
    
      axes_start_index = axes_end_index + 1;      
      
    end  

    tract_start_index = tract_end_index + 1;

  end  
  
end

  
  
