function [elems, properties] = read_fibres(filename, varargin)

%   if nargin == 2
%     num_elems = varargin{1} * 3;
%   else
%     num_elems = inf;
%   end

	file = fopen(filename,'r');	
	
	[fname, mode, machine_format] = fopen(file);

	
	if (file == -1) 
		error([ 'Could not open file ' filename '!' ]);
  end

  
	line = fgetl(file);

	if (length(line) < 13 || ~strcmp(line(1:13),'mrtrix tracks'))
		error(['File, ' filename ' was not a valid mrtrix tracks file']);
	end
	
	offset = 0;
	datatype = [];
  
  while 1
    line = fgetl(file);
    if ~ischar(line), break, end;
    line = strtrim(line);
    if strcmp(line, 'END'), break, end;
    delim = strfind (line,':');
    
    if isempty(delim)
      disp (['invalid line in header: ''' line ''' - ignored']);
    else
      
      key = lower(strtrim(line(1:delim(1)-1)));
      value = strtrim(line(delim(1)+1:end));

      if strcmp(key, 'datatype')
        
        if strcmp(value, 'Float32LE')
          machine_format = 'ieee-le';
          
        elseif strcmp(value, 'Float32BE')
          machine_format = 'ieee-be';
          
        elseif strcmp(value, 'Float64LE')
          machine_format = 'ieee-le.l64';
          
        elseif strcmp(value, 'Float64BE')
          machine_format = 'ieee-be.l64';        
          
        else
          error(['Unrecognised format ''' value '''.']);
        end

      elseif strcmp(key, 'file')
        offset = str2double(value(3:end));
      else 
        key = strrep(key, ' ', '_');
        eval(['properties.' key ' = value;']);
      end
      
    end
  end	
	
	if (~offset)
		error('"file" property was not found after 100 lines');
  end

  % Read elements
	fseek(file, offset, 'bof');
  
  if (strcmp(machine_format,'ieee-le') || strcmp(machine_format,'ieee-be'))
    elems = fread(file, inf, 'float32=>double', machine_format);
  else
    elems = fread(file, inf, 'double=>double', machine_format);
  end
%   elems = elems(1:end-3);
  
  num_triples = size(elems,1)/3;
  
  if num_triples ~= round(num_triples)
    error (['Number of elements loaded from mrtrix file data was not divisible by 3 (' num2str(size(elems,1)) ')']);
  end
  
  elems = reshape(elems,3,num_triples)';
  
  
  % Find end of file marker.
  is_pos_inf = isinf(elems) .* (elems > 0);
  
  is_pos_inf_row = is_pos_inf(:,1) .* is_pos_inf(:,2) .* is_pos_inf(:,3);
  
  end_of_data = find(is_pos_inf_row);
  
  if size(end_of_data,1) > 2
    error ([ num2str(size(end_of_data,1)) ' end_of_data markers found in file (max allowed 2). (Very) Probably infinity valued parameters have been written to fibres file.']);
    
  elseif size(end_of_data,1) == 2
    if end_of_data(2) ~= size(elems,1)
      error ('Second end_of_data marker found not at end of file, possibly corrupt file');
    end
    
    end_of_data = end_of_data(1);
  end
    
  elems = elems(1:(end_of_data-1),:);
  
end