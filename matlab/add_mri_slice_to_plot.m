function add_mri_slice_to_plot(image_location, slice_pos)

    img = load_image(image_location);
    
    data = permute(img.data(:, :, :, 1), [2, 1, 3]);
    offset = img.transform(1:3,4)';
    img_limit = offset + img.vox(1:3) .* (img.dim(1:3) - 1);
    [X,Y,Z]=meshgrid(offset(1):img.vox(1):img_limit(1),...
                     offset(2):img.vox(2):img_limit(2),...
                     offset(3):img.vox(3):img_limit(3));
    
    if any(slice_pos > img_limit) || any(slice_pos < offset)
        error(['Slice positions (' num2str(slice_pos) ') out of bounds of the image (' num2str(offset) ' -> ' num2str(img_limit) ').'])
    end
    
    hold on;
    colormap(gray)
    h = slice(X,Y,Z,data,slice_pos(1),slice_pos(2),slice_pos(3));
    set(h,'edgecolor', 'none');
    hold off;
    
    
end

