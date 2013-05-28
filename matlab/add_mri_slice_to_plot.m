function add_mri_slice_to_plot(image_location, x_slice, y_slice, z_slice)

    img = load_image(image_location);

    data = permute(img.data(:, :, :, 1), [2, 1, 3]);
    offset = img.transform(1:3,4)';
    img_limit = offset + img.vox(1:3) .* (img.dim(1:3) - 1);
    
    if any(x_slice > img_limit(1)) || any(x_slice < offset(1))
     error(['Slice positions (' num2str(x_slice) ') out of bounds of the image (' num2str(offset(1)) ' -> ' num2str(img_limit(1)) ').'])
    end

    if any(y_slice > img_limit(2)) || any(y_slice < offset(2))
     error(['Slice positions (' num2str(y_slice) ') out of bounds of the image (' num2str(offset(2)) ' -> ' num2str(img_limit(2)) ').'])
    end

    if any(z_slice > img_limit(3)) || any(z_slice < offset(3))
     error(['Slice positions (' num2str(z_slice) ') out of bounds of the image (' num2str(offset(3)) ' -> ' num2str(img_limit(3)) ').'])
    end

    [X,Y,Z]=meshgrid(offset(1):img.vox(1):img_limit(1),...
                     offset(2):img.vox(2):img_limit(2),...
                     offset(3):img.vox(3):img_limit(3)); 
    
    hold on;
    colormap(gray)
    h = slice(X,Y,Z,data,x_slice,y_slice,z_slice);
    set(h,'edgecolor', 'none');
    hold off;
    
    
end

