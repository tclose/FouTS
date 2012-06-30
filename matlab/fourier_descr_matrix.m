function psi = fourier_descr_matrix(num_samples, degree, include_endpoints) 


  if include_endpoints
    k = [0:1:(num_samples-1)]'./(num_samples-1);
  else
    k = [1:num_samples]'./(num_samples+1);    
  end

  if (degree >= 0)
      psi = ones(num_samples,1);
  end

  for (d = [1:(degree-1)])

    psi = [psi, sqrt(2) * cos(pi * k * d)];

  end

end