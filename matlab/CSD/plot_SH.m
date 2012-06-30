function plot_SH (SH, scheme, varargin)

% function plot_SH (SH, scheme)
%
% plot surface given by SH coefficients supplied in 'SH' along 
% the set of directions in 'scheme'

if (nargin < 3)
	plot_amp (SH2amp (SH, scheme), scheme);
else
	plot_amp (SH2amp (SH, scheme), scheme, [ 1 1 0 ], 1, varargin{1});
end