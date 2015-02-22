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
            for i=1:length(unique_indices)-1
                colour_indices((unique_indices(i)+1):(unique_indices(i+1))) = i;
            end
            colour_indices(last_required_index) = length(unique_indices);
        else
            colour_indices = 1:last_required_index;
        end
    else
        if length(colour_indices) < max(bundle_indices)
            error(['Length of colour indices (', num2str(length(colour_indices)),...
                   ') must be greater than or equal to the largest bundle index (',...
                  num2str(max(bundle_indices)), ')'])
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
            error(['The number of supplied colours (', num2str(length_supplied_colours),...
                   ') was less than the number required (', num2str(num_colours_required), ')']);
        end
        colours = rand([num_colours_required, 3]);
        colours_of_bundles = colours;
    end

end