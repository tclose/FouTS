function histo = lhist(gamma, mean, num_samples, range, bins)

samples = tan(pi*(rand(num_samples,1)-0.5)) * gamma + mean;

samples = samples(find(abs(samples-mean) < range));

histo = hist(samples, [-range:2*range/bins:range]);

