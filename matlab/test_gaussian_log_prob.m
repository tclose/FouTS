function log_probs = test_gaussian_log_prob(points, relative_scale)

  num_points = size(points,1);

  log_probs = zeros(num_points,1);

  for point_i = 1:num_points
    
    log_probs(point_i) = points(point_i,1).^2 + relative_scale * points(point_i,2).^2;
    
  end
  
end