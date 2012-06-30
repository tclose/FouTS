function [ F_SH, num_it ] = csdeconv(R_RH, S_SH, scheme, lambda, tau)

% function F_SH = SR_csdeconv(R_RH, S_SH, lambda, tau, scheme)
%
% Constrained spherical deconvolution:
%
% deconvolves the axially symmetric response function 'R_RH'
% (in RH coefficients) from the surface 'S_SH' (in SH 
% coefficients) to give the surface 'F_SH' (in SH coefficients),
% by constraining the corresponding amplitudes of 'F_SH' 
% calculated along the directions given by 'scheme') to be 
% non-negative. The optional parameters 'lambda' and 'tau' 
% correspond to the regularisation weight (1 by default) and the
% threshold on the FOD amplitude used to identify negative lobes
% (0.1 by default).



if ~exist('lambda'), lambda = 1; end
if ~exist('tau'), tau = 0.1; end


% guess appropriate value of lmax (assuming no super-resolution):
lmax = lmax_for_SH (S_SH);
if scheme.lmax < lmax
  lmax = scheme.lmax; 
  S_SH(nSH_for_lmax(lmax)+1:end,:) = [];
end


% build forward spherical convolution matrix up to lmax:
fconv = [];
for l = 0:2:lmax
  fconv = [ fconv; R_RH(l/2+1)*ones(2*l+1,1) ];
end
fconv = [ diag(fconv) zeros(size(fconv,1),size(scheme.sh,2)-size(fconv,1)) ];


% generate initial FOD estimate, truncated at lmax = 4:
F_SH = fconv\S_SH;
F_SH (nSH_for_lmax(4)+1:end,1) = 0;


% set threshold on FOD amplitude used to identify 'negative' values:
threshold = tau*mean (scheme.sh*F_SH);


% scale lambda to account for differences in the number of 
% SH coefficients and number of mapped directions;
lambda = lambda * size(fconv,1)*R_RH(1)/size(scheme.sh,1);


% main iteration loop:
k = [];
for num_it = 1:50
  A = scheme.sh*F_SH;
  k2 = find (A < threshold);
  if size(k2,1) + size (fconv,1) < size(scheme.sh,2)
    disp ('too few negative directions identified - failed to converge');
    return; 
  end
  if size(k) == size(k2), if k == k2, return; end; end
  k = k2;
  M = [ fconv; lambda*scheme.sh(k,:) ];
  S = [ S_SH; zeros(size(k,1),1) ];
  F_SH = M\S;
end

disp ('maximum number of iterations exceeded - failed to converge');
