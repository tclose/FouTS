function d_psi = d_fourier_descr_matrix(num_samples, degree, include_endpoints) 


  if include_endpoints
    k = [0:1:(num_samples-1)]'./(num_samples-1);
  else
    k = [1:num_samples]'./(num_samples+1);
  end

  if (degree >= 0)
      d_psi = zeros(num_samples,1);
  end

  for (d = [1:(degree-1)])

    d_psi = [d_psi, - d * pi * sqrt(2) * sin(pi * k * d)];

  end

end