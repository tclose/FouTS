function SH = gen_delta(el, az, lmax)

% function SH = gen_delta(el, az, lmax)
%
% Generate the SH coefficients for a delta function pointing
% along [ el az ] up to harmonic order 'lmax'.


SH = [];

for l = 0:2:lmax
  SH = [ SH; eval_SH(l, el, az) ];
end
