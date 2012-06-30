function [dim] = get_observed_dims(properties)

  if isfield(properties, 'exp_num_segments')
    obs_image_location = properties.obs_image;
  else
    disp('Warning!! Did not observed image location field of using number of voxels default (3).');
    num_length_sections = 3;
  end
  
  f = fopen (filename, 'r');
  if (f<1) 
    disp (['Warning!! could not open observed image file ''' filename ''' to determine number of voxels']);
    return
  end
  L = fgetl(f);
  if ~strncmp(L, 'mrtrix image', 12)
    fclose(f);
    disp ([filename ' is not in MRtrix format']);
    return
  end
  
  dim = [];
  
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
      dim = str2num(char(split_strings (value, ',')))';
      break;
    end
  end

  fclose(f);

  
  if isempty(dim)
    
  
    
end
