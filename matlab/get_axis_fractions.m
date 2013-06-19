function [ax_fractions] = get_axis_fractions(num_strands)

    dirname =fileparts(mfilename('fullpath'));
    precalc_samples_location = [dirname '/precalc_axis_fractions.mat'];
    if exist(precalc_samples_location, 'file')
       precalc_samples = load(precalc_samples_location);
    else
       precalc_samples = cell(0);
    end

    if num_strands > length(precalc_samples)
        for n=length(precalc_samples):num_strands
            precalc_samples{n} = sample_evenly_on_disc(n);
        end
        save(precalc_samples_location, precalc_samples)
    end
    ax_fractions = precalc_samples{num_strands};

end

function x = sample_evenly_on_disc(N)
    % Create randomly distributed points within a circle
    r = rand(N,1);
    theta = 2 * pi * rand(N,1);
    x0 = [cos(theta) .* r, sin(theta) .* r];
    % Perform the optimisation
    options = optimset('Algorithm', 'active-set', 'MaxFunEvals', 1e5,...
                        'MaxIter', 1e3, 'TolCon', 1e-3);
    x = fmincon(@fmin, x0, [], [], [], [], [], [], @fcon, options);
end

function y = fmin(x)
    X1 = repmat(x(:,1), [1, size(x,1)]);
    X2 = repmat(x(:,2), [1, size(x,1)]);
    d = (X1 - X1').^2 + (X2 - X2').^2;
    Y = triu(1./d.^4, 1);
    y = sum(Y(:));
end

function [y, w] = fcon(x)
    y = sum(x.^2,2) - 1;
    w = zeros(size(x));
end