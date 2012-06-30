N = 10000;
sigma_incr = 0.1;
sigma_N = 5;
output_path = '/data/home/tclose/Documents/Tractography/florey_student_talk/images/normal';
resolution = 300;
container_size = 0.35;


close all;
fig_num = 1;

x_limits = [-3.5*sigma_incr*sigma_N, 3.5*sigma_incr*sigma_N];

y_limits = [0, 1/(sigma_incr*sqrt(2*pi))];

background = [0 1 0];



%set(gca, 'visible', 'off');

% xlabel('');
% 
% colormap('gray');
% caxis([-2 -1]);
background = [0 1 0];
whitebg(background);
set(gcf, 'Color', background);
% set(gca,'XColor', [1 1 1]);
% set(gca,'YColor', [1 1 1]);
% set(gca,'ZColor', [1 1 1]);
% 
% set(gca, 'FontSize', 18);

set(gcf,'PaperUnits','centimeters');

set(gcf, 'PaperSize', [8 8]);
set(gcf, 'PaperPosition', [0 0 8 8]); 

%set(gca, 'visible', 'off');
set(gca, 'Outerposition', [0 0 1 1]);
set(gca, 'Position', [0.1 0.15 .8 .8]);

set(gcf,'InvertHardCopy','off'); 

count = 1;

range = 3.5*sigma_incr*sigma_N;
incr = range/N;

k = [-range:incr:range];

% k = [-k(end:-1:1), 0, k];

k_index = [0:1:N-1];
  
container_size_N = round(container_size * (N/range));

odd_index = k_index(find(mod(k_index ./ container_size_N, 2) >= 1))+1;

even_index = k_index(find(mod(k_index ./ container_size_N, 2) < 1))+1;

k_index(even_index) = mod(k_index(even_index), container_size_N);

k_index(odd_index) = container_size_N - mod(k_index(odd_index), container_size_N)-1;

k_index = k_index + 1;

k_index = [-k_index(end:-1:1),0,k_index] + container_size_N+1;

folded_k = [-container_size:incr:container_size];


for sigma = sigma_incr:sigma_incr:(sigma_incr*sigma_N)

  gauss_kernel = 1/(sigma*sqrt(2*pi)) * exp(- (k.^2) / (2 *sigma^2));
  
  folded_gauss_kernel = zeros(1,2*container_size_N+1);
  
  for i = 1:size(k_index,2)
    folded_gauss_kernel(k_index(i)) = folded_gauss_kernel(k_index(i)) + gauss_kernel(i);
  end
  
  
%   figure(fig_num); fig_num = fig_num + 1;
%   plot(k, gauss_kernel, 'color', [0.7383,0.6523,0.9570], 'LineWidth', 2);
  plot(folded_k, folded_gauss_kernel, 'color', [0.7383,0.6523,0.9570], 'LineWidth', 2);  
  
   xlim(x_limits);
   ylim(y_limits);
  
  
  axis('off');
  
  path = [output_path '-' num2str(count)];
  
  print(gcf, '-dpng', ['-r' num2str(resolution) ' '], path); 

  cdata = imread([path '.png']);
  % 
  imwrite(cdata, [path '.png'], 'png', 'BitDepth', 16, 'transparency', background) 

  count = count + 1;

end


return;


init_img = ones(N,N,3) * 0.5;

figure(fig_num); fig_num = fig_num + 1;
image(init_img);
pause;

b2w = -.5:(1/(N-1)):.5;
grad_img = repmat(b2w, N, 1);
grad_img = cat(3, grad_img, grad_img, grad_img);

img = init_img + grad_img;

figure(fig_num); fig_num = fig_num + 1;
image(img);
pause;

img = init_img - grad_img;

figure(fig_num); fig_num = fig_num + 1;
image(img);
pause;








for sigma = sigma_incr:sigma_incr:(sigma_incr*sigma_N)

  Sigma = sigma * N;
  k = -3*Sigma:1:(3*Sigma+1);
  
  ext_b2w = (-.5 - 3.1 * sigma):(1/(N-1)):(.5 + 3.1 * sigma);  
  
  gauss_kernel = 1/(Sigma*sqrt(2*pi)) * exp(- (k.^2) / (2 * Sigma^2));

  diff_b2w = conv(gauss_kernel, ext_b2w);
  
  trim_N = (size(diff_b2w,2) - N)/2;
  
  diff_b2w = diff_b2w(floor(trim_N):(end-ceil(trim_N)-1));

  diff_img = repmat(diff_b2w, N, 1);
  diff_img = cat(3, diff_img, diff_img, diff_img);
  
  size(diff_img)
  size(init_img)
  img = init_img + diff_img;
  
  img(find(img > 1.0)) = 1.0;
  img(find(img < 0)) = 0;  
  
  
  figure(fig_num); fig_num = fig_num + 1;
  image(img);
  pause;
  
end

img = diff_img - grad_img;

figure(fig_num); fig_num = fig_num + 1;
image(diff_img);