function [colours, colour_indices] = set_bundle_colours(bundle_indices, colours, colour_indices, compact_colours, highlight_axes)

    global colours_of_bundles;
  
    if ~exist('compact_colours', 'var')
        compact_colours = 0;
    end
    
    if ~exist('colours', 'var')
        colours = [];
    end
       
    if isempty(colours)
        colours = colours_of_bundles;
        length_supplied_colours = 0;
    else
        length_supplied_colours = size(colours, 1);
    end
    
    if ~exist('colour_indices', 'var') || isempty(colour_indices)
        last_required_index = max(bundle_indices) + 1;
        if compact_colours
            % Create a map from bundle index to a compact set of colours
            % for the bundle indices actually used.
            colour_indices = zeros(last_required_index, 1);
            unique_indices = sort(unique(bundle_indices));
            compact_index = 1;
            non_compact_index = length(unique_indices)+1;
            for bundle_i = 0:1:(last_required_index-1)
                if ismember(bundle_i, unique_indices)
                    colour_indices(bundle_i+1) = compact_index;
                    compact_index = compact_index + 1;
                else
                    colour_indices(bundle_i+1) = non_compact_index;
                    non_compact_index = non_compact_index + 1;
                end
            end
        else
            colour_indices = 1:last_required_index;
        end
    else
        if length(colour_indices) < max(bundle_indices)
            colour_indices = [colour_indices;
                ((length(colour_indices)+1):(max(bundle_indices)+1))'];
        end
        if compact_colours ~= 0
            error(['compact_colours option is not relevant when colour_indices ',...
                   'are supplied explicitly']);
        end
    end
    
    if ~exist('highlight_axes', 'var')
        highlight_axes = 0;
    end

    if highlight_axes
        num_colours_required = max(colour_indices) * 4 ;  
    else
        num_colours_required = max(colour_indices);
    end

    if size(colours, 1) < num_colours_required
        if length_supplied_colours
            warning(['The number of supplied colours (', num2str(length_supplied_colours),...
                     ') was less than the number required (', num2str(num_colours_required), ')']);
        end
        rand_colours = rand([num_colours_required - size(colours, 1), 3]);
        colours = [colours; rand_colours];
        colours_of_bundles = colours;
    end

end