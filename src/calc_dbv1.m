
%% calc_dbv1: Find closest kernel on beta manifold in L2 distance
%
% This is also used as an initial condition to the KL-div minimization
%
function [dbv1, kernc] = calc_dbv1(kern, b)
m = b/2-1; order = length(kern);
%construct points
x = zeros(order, order-m); %each column is one point (one history kernel)
for k = (m+1):order
    for t = 1:k
        x(t, k-m) = nchoosek(k, t)*(-1)^(t+1);
    end

    for t = (k+1):order
        x(t, k-m) = 0;
    end
end

X = x - x(:, 1); X = X(:, 2:end);

kern = kern(:); bhat = kern - x(:, 1);

%QR decomposition
[Q,R] = qr(X); v_qr = R \ (Q'*bhat);

kernc = x(:, 1) + X*v_qr;
dbv1 = sqrt(sum((kern-kernc).^2));

end
