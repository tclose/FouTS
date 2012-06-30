function log_probs = test_peak_log_prob(points, peaks, roi_radius, barrier_rate)

  num_peaks = size(peaks,1);
  num_points = size(points,1);
  
  log_probs = zeros(num_points,1);

  for peak_i = 1:num_peaks

    peak_centre = peaks(peak_i, 1:2);
    peak_width = peaks(peak_i, 3:4);
    peak_height = peaks(peak_i, 5);        
    peak_type = peaks(peak_i, 6);  

    if peak_type == 0 % Gaussian

      peak_variance = (peak_width / 3.0).^2;      

      peak_log_probs = peak_height * exp (- ((peak_centre(1) - points(:,1)).^2) ./ peak_variance(1)) .* exp (- ((peak_centre(2) - points(:,2)).^2) ./ peak_variance(2));

    elseif peak_type == 1 % Pyramid

      peak_steepness(1) = peak_height / peak_width(1);
      peak_steepness(2) = peak_height / peak_width(2);      

      peak_log_probs = peak_height - abs(peak_centre(1) - points(:,1)) .* peak_steepness(1) - abs(peak_centre(2) - points(:,2)) .* peak_steepness(2);

      peak_log_probs(find(peak_log_probs < 0)) = 0;

    else
      error(['Unrecognised peak type, ' num2str(peak_type) '.']);
    end

    points_norm = sqrt(points(:,1).^2 + points(:,2).^2);
    
    outside_roi = find(points_norm > roi_radius);
    
    peak_log_probs(outside_roi) = peak_log_probs(outside_roi) - ((points_norm(outside_roi) - roi_radius) * barrier_rate).^2;

    log_probs = log_probs + peak_log_probs;

  end
  
end