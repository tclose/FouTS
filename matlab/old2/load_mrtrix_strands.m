function strands = load_mrtrix_strands(filename, varargin)
% samples = f_descript(num_samples, coeff, [plot_type])
% 
%	Args:			
%
%			samples - the samples, which are plotted
%			num_samples - number of samples
%			plot_type (optional - passed directly to plot3 (e.g. 'x' plots
%												crosses instead of a line)

	if (nargin > 4)
		error(['Incorrect number of arguments ' nargin ' expecting no more than 2']);
		
	elseif (nargin == 3) 
		first_strand = 1;
		last_strand = varargin{2};
		
	elseif (nargin == 4) 
		first_strand = varargin{2};
		last_strand = varargin{3};
	else
		first_strand = 1;
		last_strand = -1;
	end
	
	if (nargin == 2)
		strand_diameter = varargin{1};
	else
		strand_diameter = 0.03;
	end
	
	dots_i = findstr(filename,'.');
	ext = filename(dots_i(end):end);
	
	convert_to_tck = strcmp(ext, '.frr');
		
	
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
	datatype = [];
		
	for (line_i = 1:20)	
		line =  fgetl(file);
		if (strcmp(line(1:4),'file'))
			offset = str2double(line(9:end));
			break
% 		elseif (strcmp(line(1:8),'datatype'))
% 			datatype = line(11:end);
 		end
	end;	
	
	if (~offset)
		error(['"offset" property was not found after 10 lines']);
	end
	
% 	if (datatype == [])
% 		error(['"datatype" property was not found after 10 lines']);
% 	end


	fseek(file, offset, 'bof');

	elems = fread(file, inf, 'float32=>double', machine_format);
	
	
	end_of_strands = find(isnan(elems));
	
	end_of_strands = end_of_strands(1:3:end) - 1;
%  	for (i = 1: length(end_of_strands))
%  		disp([ num2str(end_of_strands(i)) ', ' num2str(elems(end_of_strands(i))) ', ' num2str(elems(end_of_strands(i)+1)) ', ' num2str(elems(end_of_strands(i)+2)) ', ' num2str(elems(end_of_strands(i)+3)) ', ' num2str(elems(end_of_strands(i)+4))])  
% 	end
	
	total_number_of_strands = size(end_of_strands,1);
	
	if (last_strand == -1 | (last_strand > total_number_of_strands) )
		last_strand = size(end_of_strands,1);
	end
	
	number_of_strands = last_strand - first_strand + 1;
	
	if (first_strand == 1)
		end_of_prev_strand = -3;
	else
		end_of_prev_strand = end_of_strands(first_strand-1);
	end
	
	strands = cell(number_of_strands,4);
	
	for (strand_i = first_strand:last_strand) 
		
		output_strand_i = strand_i - first_strand + 1;
		
		strand = elems(end_of_prev_strand+4:end_of_strands(strand_i));
		
		num_points = size(strand,1) / 3;
		
		if (num_points ~= floor(num_points))
			error(['Strand ' strand_i ' does not contain a number of elements that are divisible by 3, ' num2str(size(strand,1))]);
		end
		
		strand = reshape(strand,3,num_points)';
		
		if (convert_to_tck)
			strand = fourier2tck(strand);
		end
		
		strands{output_strand_i,1} = strand;
		strands{output_strand_i,2} = strand_i-1;
		strands{output_strand_i,3} = strand_diameter;
		strands{output_strand_i,4} = strand_i-1;
		
		end_of_prev_strand = end_of_strands(strand_i);
  end
  
  strands{end+1,1} = filename;
	
end

