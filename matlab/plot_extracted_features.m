function plot_extracted_features(directory)

global v0_norm v1_norm v2_norm dot_v1_v2 dot_v0_v1 angle_v1_v2

load([ directory '/angle_v1_v2.txt']);

load([ directory '/dot_v0_v1.txt']);
load([ directory '/dot_v1_v2.txt']);
load([ directory '/v0_norm.txt']);
load([ directory '/v1_norm.txt']);
load([ directory '/v2_norm.txt']);

figure_start = gcf + 1;

figure(figure_start+1);
hist(v0_norm,1000)
title('v0 norm');

figure(figure_start+2);
hist(v1_norm,1000)
title('v1 norm');

figure(figure_start+3);
hist(v2_norm,1000)
title('v2 norm');

figure(figure_start+4);
hist(dot_v1_v2,1000)
title('dot v1 v2');

figure(figure_start+5);
hist(dot_v0_v1,1000)
title('dot v0 v1');

% figure(6);
% hist(angle_v1_v2,1000)
% title('angle v1 v2');

figure(figure_start+7);
hist2d(dot_v1_v2, v2_norm,100);
title('dot v1 v2 - v2 norm');

figure(figure_start+8);
hist2d(dot_v0_v1, v0_norm,100);
title('dot v0 v1 - v0 norm');


end