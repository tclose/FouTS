function [dims, true_location, voxel_lengths, offset] = get_observed_properties(properties)

  true_location = [];

  if isfield(properties, 'obs_image')
    obs_image_location = properties.obs_image;
  else
    disp('Warning!! Did not find observed image location field of using  image dimensions default (3x3x3).');
    dims = [0 0 0];
    voxel_lengths = [0 0 0 0];
    offset = [0 0 0];
    return;
  end
  
  f = fopen (obs_image_location, 'r');
  if (f<1)
    if strcmp(obs_image_location(1:6),'/work/')
        index = strfind(obs_image_location, '/params/');
        f = fopen(strcat(['/home/tclose/fouts/', obs_image_location(index:end)]));
    end
    if (f<1)
        disp (['Warning!! could not open observed image file ''' obs_image_location ''' to determine image dimensions']);
        dims = [3,3,3];
        return
    end
  end
  L = fgetl(f);
  if ~strncmp(L, 'mrtrix image', 12)
    fclose(f);
    disp (['Warning!! observed image file ''' obs_image_location ''' is not in MRtrix format, could not determine image dimensions.']);
    dims = [3 3 3];
    return
  end
  
  dims = [];
  voxel_lengths = [];
  offset = {};
  
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
      elseif strcmp(key, 'vox')
        voxel_lengths = str2num(char(split_strings (value, ',')))';        
      elseif strcmp(key, 'transform')
        row = str2num(char(split_strings (value, ',')))'; %#ok<*ST2NM>
        offset{end+1} = row(4);  %#ok<AGROW>
      elseif strcmp(key, 'state_location')
        true_location = char(split_strings (value, ','));
      end
    end
  end
  
  offset = cell2mat(offset);
  
  fclose(f);
 
%   if isempty(dims)
%     disp('Warning!! Observed image location did not specifiy image dimensions.');
%     dims = [3 3 3];
%   end
%   
%   if isempty(true_location)
%     disp('Warning!! Observed image location did not specifiy inital location.');
%   end

  
    
end

  

function S = split_strings (V, delim)
  S = {};
  while size(V,2) > 0
    [R, V] = strtok(V,delim);
    S{end+1} = R;
  end
end