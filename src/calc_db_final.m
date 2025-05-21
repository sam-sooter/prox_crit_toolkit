function [db, sddb, kern, sigma, H, kernc, exit_status] = calc_db_final(x, order, deltaT, b, with_err_bars, with_QC, with_parallel, fit_method)

%% INPUTS

% x: the observed time series. If the time series is segmented, store the
% segments in elements of a cell array e.g. x = {[1, 2, 3], [10, 5, 0]}
% means the time series consists of two separate segments, [1, 2, 3] and
% [10, 5, 0]

% order (int): the AR model order

% deltaT (double): the time between consecutive points in the time series, in units
% of seconds

% b (int): the type of criticality to calculate proximity to (b = 2, 4, 6, ...,
% 2*order)

% with_err_bars (bool): if true, calculate error bars for db

% with_QC (bool): if true, make quality control plots

% with_parallel (bool): if true, use parallel computing (only speeds up
% calculation if with_err_bars = true

% fit_method (string): 'YuleWalker' to use Yule-Walker method for AR model
% fitting, 'MaxLikelihood' to use maximum likelihood AR fitting

%% OUTPUTS

% db (double): proximity to criticality of type b in units of bits/s

% sddb (double): error bar for db (nan unless with_err_bars = true and 
% fit_method = 'MaxLikelihood')

% kern (array): history kernel of best-fit AR model

% sigma (double): noise standard deviation of best-fit AR model

% H (2d array): hessian of log-likelihood function (nan unless fit_method = 
% 'MaxLikelihood')

% kernc (array): history kernel of closest (in KL div. rate) AR model in
% the b critical manifold

% exit_status (int): 0 if all calculations were successful, 1 otherwise


exit_status = 0;

if ~iscell(x); x = {x}; end

% z-score
tot = [];
for i = 1:length(x); tot = [tot x{i}(:)']; end
for i = 1:length(x); x{i} = (x{i}-mean(tot))/std(tot); end

if strcmp(fit_method, 'MaxLikelihood')
    % max likelihood fit
    [kern, sigma, H] = mle_ar(x, order);
elseif strcmp(fit_method, 'YuleWalker')
    % yule walker fit
    [kern, sigma] = yw_ar(x, order);
    H = nan;
else
    exit_status = 1;
    db = nan; sddb = nan; kern = nan; sigma = nan; H = nan; kernc = nan;
    disp("Error: invalid fit_method. Change to 'MaxLikelihood' or 'YuleWalker'");
    return;
end

% check whether best-fit model is explosive
if explosive(kern)
    exit_status = 1;
    db = nan; sddb = nan; kernc = nan; disp('Error: estimated model is explosive.')
    return;
end

% calculate db
[db, cm, exit] = minimize_klr_b(kern, b);
if exit == 1; exit_status = exit; end

% closest critical AR model
kernc = cm.kernc;

% change units to bits/s
db = log2(exp(1))/deltaT*db;

% quality control plot (check optimization)
if with_QC
    QC_plot(kern, kernc, db, b, deltaT);
end

% calculate error bars for db
if with_err_bars
    % fit_method must be 'MaxLikelihood' to calculate error bars
    if strcmp(fit_method, 'YuleWalker')
        exit_status = 1; sddb = nan;
        disp('Error: error bar calculation not supported with YuleWalker fitting. Switch to MaxLikelihood if you want error bars.')
        return;
    end
    % gradient of db with respect to kern
    graddb = zeros(order, 1);

    % step size for numerical calculation of gradient
    STEP = 0.001;
    check = true;
    while check
        for j = 1:order
            % if STEP will remove kern from stable region, reduce STEP
            if explosive(kern(:)+STEP*[zeros(1, j-1) 1 zeros(1, order-j)]')
                STEP = STEP/10;
                break;
            end
        end
        check = false;
    end

    if with_parallel
        exit_temp = zeros(1, order);
        parfor j = 1:order
            % calculate db at perturbed kernel kern + STEP*[0 ... 1 ... 0]
            [newdb, ~, exit_temp(j)] = minimize_klr_b(kern(:)+STEP*[zeros(1, j-1) 1 zeros(1, order-j)]', b);
       
            % numerical estimate of gradient
            newdb = log2(exp(1))/deltaT*newdb;
            graddb(j) = (newdb-db)/STEP;
        end
        if sum(exit_temp)>0; exit_status = 1; end
    else
        for j = 1:order
            % calculate db at perturbed kernel kern + STEP*[0 ... 1 ... 0]
            [newdb, ~, exit] = minimize_klr_b(kern(:)+STEP*[zeros(1, j-1) 1 zeros(1, order-j)]', b);
            if exit == 1; exit_status = 1; end
            % numerical estimate of gradient
            newdb = log2(exp(1))/deltaT*newdb;
            graddb(j) = (newdb-db)/STEP;
        end
    end

    % error bars for db
    sddb = sqrt(-graddb'*inv(H)*graddb);
else
    sddb = nan;
end
end

% checks whether the AR model with history kernel kern is explosive
function expl = explosive(kern)
ERROR = 1e-2;
expl = min(abs(roots([fliplr(kern(:)') -1])))<1-ERROR;
end

% fourier transforms kern
function kernw = get_kernw(kern, w)
kernw = zeros(size(w));
for t = 1:length(kern)
    kernw = kernw + kern(t)*exp(-1i*w*t);
end
end

% calculate the value of sigmac that minimizes KL div. rate for a given
% kern and kernc
function sigmac = get_sigmac(kern, kernc)
g = @(w) abs(1-get_kernw(kern, w)).^2;
gc = @(w) abs(1-get_kernw(kernc, w)).^2;
fun = @(w) 1/(2*pi)*gc(w)./g(w);
sigmac = sqrt(integral(fun, -pi, pi, 'AbsTol', 1e-10));
end

% calculate KL div. rate between 
function klr = get_klr_b(kern, kernc, b)
kernc = beta_helper(kernc, b); 
sigmac = get_sigmac(kern, kernc);
g = @(w) abs(1-get_kernw(kern, w)).^2;
gc = @(w) abs(1-get_kernw(kernc, w)).^2;
fun2 = @(w) 1/(4*pi)*(log(sigmac^2*g(w)./gc(w))+gc(w)./(sigmac^2*g(w))-1);
xx = linspace(-pi, pi, 1000);
klr = trapz(xx, fun2(xx));
%klr = integral(fun2, -pi, pi, 'AbsTol', 1e-10);
end

% constraint
function [c,ceq] = con(kernc, b)
    kernc = beta_helper(kernc, b); kernc = kernc(:)';
    c = 1-min(abs(roots([fliplr(kernc(:)') -1])));
    ceq = [];
end

% minimize KL div rate over the critical manifold
function [db, closest_model, exit] = minimize_klr_b(kern, b)
exit = 0;
closest_model = struct;
if length(kern)<b/2
    exit = 1;
    db = nan; disp('WARNING: invalid model order');
elseif length(kern)==b/2
    db = get_klr_b(kern, [], b);
    temp = beta_helper([], b);
    closest_model.sigmac = get_sigmac(kern, temp);
    closest_model.kernc = temp;
else
    objective = @(kernc) get_klr_b(kern, kernc, b);
    
    % Starting point
    [~, cme] = calc_dbv1(kern, b);
    kernc0 = cme(1:(length(kern)-b/2));
    
   
    % Call fmincon
    options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp',  'ConstraintTolerance', 0.01, 'MaxFunctionEvaluations', 10000, 'StepTolerance', 1e-15);

    nonlcon = @(kernc) con(kernc, b);
    [x, db] = fmincon(objective, kernc0, [], [], [], [], [], [], nonlcon, options);

    closest_model.kernc = beta_helper(x, b);
    closest_model.sigmac = get_sigmac(kern, closest_model.kernc);

    if explosive(closest_model.kernc)
        exit = 1; db = nan;
        disp('ERROR: drifted off of critical manifold.')
   end
end
end



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


function [phic] = beta_helper(phic, beta)
phic = phic(:);
n = length(phic)+beta/2; %model order
A = zeros(beta/2);
b = zeros(beta/2, 1);
for i = 1:size(A, 1)
    A(i, :) = ((n-size(A, 2)+1):n).^(i-1);
    if i == 1
        b(i) = 1-sum(phic);
    else
        b(i) = -((1:length(phic)).^(i-1))*phic;
    end
end

phic = [phic ; A\b];
end



% AR model fitting
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

kern = inv(Gamma)*gamma(:);
sigma = sqrt(1-gamma(:)'*kern(:));

end

function QC_plot(kern, kernc, db, b, deltaT)
figure; dim = length(kern)-b/2;
for i = 1:dim
    subplot(ceil(sqrt(dim)), ceil(sqrt(dim)), i);
    hold on;
    xx = linspace(kernc(i)-0.2, kernc(i)+0.2, 100);
    klr = zeros(size(xx));
    expl = zeros(size(xx));
    for j = 1:length(xx)
        kernc_temp = kernc(1:dim); kernc_temp(i) = xx(j);
        kernc_temp = beta_helper(kernc_temp, b);
        klr(j) = log2(exp(1))/deltaT*get_klr_b(kern, kernc_temp(1:dim), b);
        expl(j) = explosive(kernc_temp);
    end
    plot(xx(expl==1), klr(expl==1), '.r', 'MarkerSize', 18);
    plot(xx(expl==0), klr(expl==0), '.c', 'MarkerSize', 18);

    xline(kernc(i), 'k', 'LineWidth', 2)
    yline(db, 'k', 'LineWidth', 2)
    xlabel(['$(\phi_c)_' num2str(i) '$'], 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('KL rate (bits/s)', 'FontSize', 12);
    set(gca, 'LineWidth', 1); box on;
end
end