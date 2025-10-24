
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AR Model fitting
%% 
%% mle_ar
% Use maximum likelihood to fit an AR model to x
function [kern, sigma, H] = mle_ar(x, order)
if ~iscell(x); x = {x}; end
T = 0; for i = 1:length(x); T = T + length(x{i}); end
X = zeros(T, order); % rows are different samples
Y = zeros(T, 1);
counter = 1;
for i = 1:length(x)
    for j = 1:length(x{i})
        Y(counter) = x{i}(j);
        hist = x{i}(max([1, j-order]):(j-1));
        X(counter, :) = [fliplr(hist(:)') zeros(1, order-length(hist))];
        counter = counter + 1;
    end
end

kern = inv(X'*X)*X'*Y;
sigma = sqrt(mean((Y-X*kern).^2));

H = zeros(order);
A = sum((Y-X*kern).^2);
B = zeros(order, 1);
for j = 1:order
    B(j) = sum((Y-X*kern).*X(:, j));
end

for j = 1:order
    for k = 1:order
        H(j, k) = T*(2*B(j)*B(k)/A^2-1/A*sum(X(:, j).*X(:, k)));
    end
end
end