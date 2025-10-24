

%% yw_ar
% Use the Yule-Walker method to fit AR model to data
%
function [kern, sigma] = yw_ar(x, order)
% estimate autocovariance up to lag order
% x is already z-scored (variance = 1, mean = 0)
gamma = zeros(1, order);
for lag = 1:order
    cntr = 0;
    for j = 1:length(x)
        for k = 1:(length(x{j})-lag)
            gamma(lag) = gamma(lag) + x{j}(k)*x{j}(k+lag);
            cntr = cntr + 1;
        end
    end
    gamma(lag) = gamma(lag)/cntr;
end

Gamma = ones(order);
for i = 1:order
    for j = (i+1):order
        Gamma(i, j) = gamma(abs(i-j));
        Gamma(j, i) = Gamma(i, j);
    end
end

kern = Gamma\gamma(:);
sigma = sqrt(1-gamma(:)'*kern(:));


end

