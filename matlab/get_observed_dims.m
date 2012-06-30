function [dims, true_location] = get_observed_properties(properties)

  if isfield(properties, 'obs_image')
    obs_image_location = properties.obs_image;
  else
    disp('Warning!! Did not find observed image location field of using  image dimensions default (3x3x3).');
    dims = [3 3 3];
    return;
  end
  
  f = fopen (obs_image_location, 'r');
  if (f<1) 
    disp (['Warning!! could not open observed image file ''' obs_image_location ''' to determine image dimensions']);
    dims = [3,3,3];
    return
  end
  L = fgetl(f);
  if ~strncmp(L, 'mrtrix image', 12)
    fclose(f);
    disp (['Warning!! observed image file ''' obs_image_location ''' is not in MRtrix format, could not determine image dimensions.']);
    dims = [3 3 3];
    return
  end
  
  dims = [];
  
  while 1
    L = fgetl(f);
    if ~ischar(L), break, end;
    L = strtrim(L);
    if strcmp(L, 'END'), break, end;
    d = strfind (L,':');
    if ~isempty(d)
      key = lower(strtrim(L(1:d(1)-1)));
      value = strtrim(L(d(1)+1:end));
      if strcmp(key, 'dim')
        dims = str2num(char(split_strings (value, ',')))';
        break;
      end
    end
  end
  
  fclose(f);

  if isempty(dims)
    disp('Warning!! Observed image location did not specifiy image dimensions.');
    dims = [3 3 3];
  end

  
    
end

  

function S = split_strings (V, delim)
  S = {};
  while size(V,2) > 0
    [R, V] = strtok(V,delim);
    S{end+1} = R;
  end
end