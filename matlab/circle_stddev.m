

samp_period = 0.001;

x = [-1:samp_period:1]';
y = [-1:samp_period:1]';

[X,Y] = meshgrid(x,y);

R = X.^2 + Y.^2;

ind = find(R>1);

X(ind) = [];

x = X(:);

% hist(x,100);

std(x)



