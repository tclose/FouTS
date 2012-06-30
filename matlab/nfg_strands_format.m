function nfg_strands = nfg_strands_format(strands, prop_keys, prop_values, filename, varargin)

	if (nargin > 6)
		error(['Incorrect number of arguments ' nargin ' expecting no more than 4']);
  end

  if (nargin >= 5)
    num_points = varargin{1};
  else  
    num_points = 100;
  end
  
  if (nargin == 6)
		strand_radius = varargin{2};
	else
		strand_radius = [];
  end

  already_track = strcmp(file_extension(filename),'tck');
  
  strand_radius_key_i = strmatch('track_radius', prop_keys, 'exact');
  bundle_index_key_i = strmatch('bundle_index', prop_keys, 'exact');

  
  num_strands = size(strands,1);
  
  nfg_strands = cell(num_strands, 4);
  
  strand_indices = [1:num_strands]';
  
  if ~isempty(strand_radius)
    
    strand_radii = ones(num_strands,1) * strand_radius;
    
  else
      
    if isempty(strand_radius_key_i)      
      error(['No strand radii loaded from extended properties, so strand radius must be provided']);
    end
    for strand_i = 1:num_strands
      strand_radii(strand_i) = str2num(prop_values{strand_i,strand_radius_key_i});
    end
  end

  if ~isempty(bundle_index_key_i)      
    for strand_i = 1:num_strands
       bundle_indices(strand_i) = str2num(prop_values{strand_i,bundle_index_key_i});  
    end
  else
    bundle_indices = strand_indices;
  end


   for strand_i = 1:num_strands
     
     if ~already_track
       nfg_strands{strand_i,1} = fourier2tck(strands{strand_i});
     else
       nfg_strands{strand_i,1} = strands{strand_i};
     end
     nfg_strands{strand_i,2} = strand_indices(strand_i);
     nfg_strands{strand_i,3} = strand_radii(strand_i);
     nfg_strands{strand_i,4} = bundle_indices(strand_i);
   end

   nfg_strands{end+1,1} = strrep(filename, '_', ' ');
   
end