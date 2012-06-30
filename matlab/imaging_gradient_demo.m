close all;

N = 100;
sigma = 0.05 * N;


kernel_extent = 3;

kernel_size = ceil(sigma * kernel_extent)

gauss_kernel = zeros(2 * kernel_size + 1);

for i = 1:size(gauss_kernel,1)
  
  ker_t = kernel_extent * (i - kernel_size - 1)/kernel_size;
 
  gauss_kernel(i) = 1/(sigma * sqrt(2) * pi ) * exp(- ker_t^2/(2 * sigma^2));
end


plot(gauss_kernel);

return;

initial = [0:1/(N-1):1];
final = [1:-1/(N-1):0];

% G = repmat(row, M, 1);


imagesc(initial);
colormap('gray');
figure;
imagesc(final);
colormap('gray');


