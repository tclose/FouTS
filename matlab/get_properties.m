function values = get_properties(keys, all_values, key, varargin)

  if nargin >= 4
    default_value = varargin{1};
       
  else
    default_value = [];
  end

  if nargin >= 5
    num_rows = varargin{2};
  else
    num_rows = 0;
  end

  if nargin == 6
    datatype = varargin{3};  % Defaults to double but can also be 'string'.
  else
    datatype = 'double';
  end

  key_i = strmatch(key, keys, 'exact');
  
  if ~isempty(key_i)
    
    if strcmp(datatype, 'double')
      values = str2double(all_values(:, key_i));
      
    elseif strcmp(datatype, 'string')
      values = all_values(:, key_i);
      
    else
      error(['Unrecognised dataype ''' datatype '''']);
      
    end
    
  else
    
    if size(default_value) == [1 1]; %If the default value is a scalar make a vector the size of the number of rows with each element equal to the default value 
      values = ones(num_rows,1) * default_value;
      
    elseif (size(default_value,1) == num_rows) || (size(keys,1) == 0) % Else if the default value
      values = default_value;
      
    else
      error(['Size of default value (' num2str(size(default_value,1)) ') does not match number of rows (' num2str(num_rows) '), and alternatively is not a scalar']);
      
    end
      
  end


end