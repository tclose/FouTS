function [mcmc_state_matrix, labels] = load_unzip_mcmc_state(filename)

  mcmc_state_matrix = load(filename)';
  
  num_elems = size(mcmc_state_matrix,1);
  
  labels = cell(num_elems,1);
  
  for state_i = 1:num_elems
    labels{state_i} = ['Elem ' num2str(state_i)];
  end
  
end