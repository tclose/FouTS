function [strands, bundle_indices] = tracts2strands(tracts, base_widths, num_strands, highlight_axes, oblong, tract_indices, strands_per_acs, acs, accurate_acs)

  if exist('num_strands', 'var')
      if num_strands < 0
          error(['''-num_strands'' must be greater than zero (' num2str(num_strands) ').']);  
      end
  else
      num_strands = 4;
  end

  if ~exist('highlight_axes','var')
    highlight_axes = 0;
  end

  if ~exist('oblong','var')
    oblong = 0;
  end
  
  if ~exist('tract_indices','var')
    tract_indices = [];
  end
  
  if ~exist('strands_per_acs','var')
    strands_per_acs = -1;
  end
  
  if ~exist('acs', 'var')
    acs = ones(length(tracts),1);
  end
  
  % Uses axis fractions generated from optimisations, allowing the ACS
  % to be represented by whole number increments in number of strands 
  % rather than whole number increments in the square-root of strands
  if ~exist('acurate_acs','var')
    if strands_per_acs < 0
      accurate_acs = false;
    else
      accurate_acs = true;
    end
  end
  
  if isempty(tract_indices)
    tract_indices = [0:1:size(tracts,1)-1];
  end
  
  num_tracts = size(tracts,1);
  
  strands = cell(num_tracts,1);
  bundle_indices = zeros(num_tracts,1);
  
  strand_count = 0;
  tcks = cell(num_tracts * (2 * num_strands + 1),4);
    
  for tract_i = 1:num_tracts

    base_index = tract_indices(tract_i);
    
    if highlight_axes
      base_index = base_index * 4;
    end
  
    if accurate_acs
       ax_fractions = get_axis_fractions(round(strands_per_acs * acs(tract_i)));
    else
        if strands_per_acs ~= -1
          num_strands = round(sqrt(strands_per_acs * acs));    
        end

        width_fractions = ones(size(tracts,1),1) ./ num_strands;

        ax_frac_col = (-1+width_fractions(tract_i)):(2*width_fractions(tract_i)):(1-width_fractions(tract_i));
        ax_fractions = [ax_frac_col; ax_frac_col];
    end
    
    if ~isempty(ax_fractions)
      for ax2_frac = ax_fractions(1, :)

        for ax3_frac = ax_fractions(2, :)

          if (sqrt(ax2_frac^2 + ax3_frac^2) <= 1.0) || oblong
            strand_count = strand_count + 1;

            strands{strand_count} = tracts{tract_i, 1} + (tracts{tract_i, 2} * ax2_frac + tracts{tract_i, 3} * ax3_frac) * base_widths(tract_i);

            if highlight_axes

              if ax2_frac == 0 && ax3_frac == 0
                strand_index = base_index + 1;
              elseif ax2_frac == 0
                strand_index = base_index + 2;
              elseif ax3_frac == 0
                strand_index = base_index + 3;
              else
                strand_index = base_index;
              end

            else
              strand_index = base_index;
            end

            bundle_indices(strand_count) = strand_index;

          end
        end
      end 
    end

  end  


end