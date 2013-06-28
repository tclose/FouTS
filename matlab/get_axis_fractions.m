function [ax_fractions] = get_axis_fractions(num_strands)
    if num_strands == 0
        ax_fractions = [];
    else
        dirname =fileparts(mfilename('fullpath'));
        precalc_samples_location = [dirname '/precalc_axis_fractions.mat'];
        if exist(precalc_samples_location, 'file')
           load(precalc_samples_location, 'precalc_samples');
        else
           precalc_samples = cell(0);
        end

        if num_strands > length(precalc_samples)
            for n=(length(precalc_samples)+1):num_strands
                precalc_samples{n} = sample_evenly_on_disc(n);
            end
            save(precalc_samples_location, 'precalc_samples')
        end
        ax_fractions = precalc_samples{num_strands}';
    end
end

function x = sample_evenly_on_disc(N)
    if N == 1
        x = [0, 0];
    else
        % Create randomly distributed points within a circle
        x0 = zeros(N,2);
        step = 2 / floor(sqrt(N));
        count = 0;
        for x=-1:step:1
            for y=-1:step:1
                if x^2 + y^2 < 1
                    count = count + 1;
                    x0(count,:) = [x, y];
                end
            end
        end 
        circ_step = 2*pi / (N - count);
        for theta = 0:circ_step:(2*pi - circ_step)
            count = count + 1;
            x0(count, :) = [cos(theta) sin(theta)];
        end
        if count ~= N
            error('count doesn''t equal N')
        end
        % Perform the optimisation
        options = optimset('Algorithm', 'active-set', 'MaxFunEvals', 1e5,...
                            'MaxIter', 1e3, 'TolCon', 1e-3);
        x = fmincon(@fmin, x0, [], [], [], [], [], [], @fcon, options);
        % Calculate the scale away from the edge of the tube
        x = x * (1 - 1 / ceil(sqrt(N)));
    end
end

function y = fmin(x)
    X1 = repmat(x(:,1), [1, size(x,1)]);
    X2 = repmat(x(:,2), [1, size(x,1)]);
    d = ((X1 - X1').^2 + (X2 - X2').^2).^4;
    Y = triu(1./d, 1);
    y = sum(Y(:));
%      d = d + eye(length(x)) * 2;
%      y = -max(min(d,[],2),[],1);
end

function [y, w] = fcon(x)
    y = sum(x.^2,2) - 1;
    w = zeros(size(x));
end