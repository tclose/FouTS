function s = eval_SH(l, el, az)

% syntax: s = eval_SH(l, el, az)
%
% Evaluates the lth order spherical harmonic coefficients at
% positions [ el az ].

s = ones(size(az,1),1);

if l > 0  
  s = [ sin(az*(l:-1:1)) s cos(az*(1:l)) ];
end

s = eval_ALP(l, el).*s';
