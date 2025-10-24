function [db, sddb, kern, sigma, H, kernc, exit_status] = calc_db(x, order, deltaT, b, with_err_bars, with_QC, with_parallel, fit_method)

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

% If x is a single time series: cast as a 1x1 cell array.
if ~iscell(x); x = {x}; end

% z-score
tot = [];
for i = 1:length(x); tot = [tot x{i}(:)']; end
for i = 1:length(x); x{i} = (x{i}-mean(tot))/std(tot); end

% Find the best-fit AR(n) model to "x"
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% minimize_klr_b
%  Main function: minimize KL div rate over the desired critical
%  manifold (Mb)
%    
%%  Inputs:
%   kern:  Kernel of an AR(n) model ("n" is determined by vector length
%
%   b (int, must be even): Order of the critical manifold Mb (2,4,6,...etc)
%
%%  Outputs:
%   db:  Distance from Mb
%
%   closest_model: structure that includes
%       kernc:   kernel for the closest model
%       sigmac:  noise std for the closest model
%
%   exit:   exit flag (0=successful)
%
%% Notes: Within this function, and any functions it calls, 
%  kernels on the critical manifold Mb are represented by a 
%  reduced kernel (which is just the first n-(b/2) entries
%  To go from 
% 
function [db, closest_model, exit] = minimize_klr_b(kern, b)
exit = 0;
closest_model = struct;
if length(kern)<b/2
    exit = 1;
    db = nan; disp('WARNING: invalid model order');
elseif length(kern)==b/2
    db = get_klr_b(kern, [], b,'integral');
    temp = beta_lift([], b);
    closest_model.sigmac = get_sigmac(kern, temp);
    closest_model.kernc = temp;
else
    objective = @(kernc) get_klr_b(kern, kernc, b,'integral');
    
    % Starting point: solve the linear system to minimze L2 
    %    distance to Mb
    [~, cme] = calc_dbv1(kern, b);

    % Project to reduced kernel
    kernc0 = beta_project(cme,b);
    
   
    % Call fmincon
    options = optimoptions('fmincon', 'Display', 'none', ...
        'Algorithm', 'sqp',  ...
        'ConstraintTolerance', 0.01, ...
        'MaxFunctionEvaluations', 10000, ...
        'StepTolerance', 1e-15);

    nonlcon = @(kernc) con(kernc, b);
    [x, db] = fmincon(objective, kernc0, [], [], [], [], [], [], nonlcon, options);

    closest_model.kernc = beta_lift(x, b);
    closest_model.sigmac = get_sigmac(kern, closest_model.kernc);

    if explosive(closest_model.kernc)
        exit = 1; db = nan;
        disp('ERROR: drifted off of critical manifold.')
   end
end
end


%% get_klr_b: calculate KL div. rate between two AR models,
%%   one of which lies on the b=beta manifold (Mb)
% 
%% Inputs:
%
%  kern:   (1xn) kernel for an AR(n) model
%
%  kernc:  (1 x (n-b/2)) reduced kernel for an AR(n) model on Mb
%
%  b: which critical manifold
%
%  int_method: currently 'trapz' or 'integral'
function klr = get_klr_b(kern, kernc, b, int_method)
kernc = beta_lift(kernc, b); 
sigmac = get_sigmac(kern, kernc);
g = @(w) abs(1-get_kernw(kern, w)).^2;
gc = @(w) abs(1-get_kernw(kernc, w)).^2;
fun2 = @(w) 1/(4*pi)*(log(sigmac^2*g(w)./gc(w))+gc(w)./(sigmac^2*g(w))-1);

if strcmp(int_method, 'trapz')
    xx = linspace(-pi, pi, 1000);
    klr = trapz(xx, fun2(xx));
elseif strcmp(int_method,'integral')
    klr = integral(fun2, -pi, pi, 'AbsTol', 1e-10);
else
    disp('ERROR (get_klr_b): specific integration method');
    disp('Valid options are "trapz" and "integral"'); 
end

end

%% con: Constraints for optimization problem
% 
%% Outputs:
%
%  c: Inequality constraints ("c<0"). Currently, this is that the zeros of
%  the AR polynomial must lie outside the unit disk.
%
%  ceq: Equality constraints ("ceq=0"). Currently none.
%   
function [c,ceq] = con(kernc, b)
    kernc = beta_lift(kernc, b); kernc = kernc(:)';
    c = 1-min(abs(roots([fliplr(kernc(:)') -1])));
    ceq = [];
end



%% beta_project: Given a point on the critical manifold,
%% project it into a reduced space
function [phic_reduc] = beta_project(phic, beta)
phic = phic(:);
n = length(phic);    %model order

if (beta/2>n)
    phic_reduc = [];
    disp('ERROR (beta_project): beta is too large '); 
else
    phic_reduc = phic(1:(n-beta/2));

end
end

%% beta_lift: find a point on critical manifold Mb, given 
%%     a reduced kernel
%
%% Inputs:
%
% phic_reduc ((1 x n-beta/2) array:   first n-(beta/2) kernel coefficients
%
% beta:  order of the critical manifold
%
%% Outputs: 
% phic:  ((1 x n) array:   full length n kernel.
%   Uniquely determined by the input kernel, by applying linear 
%   constraints outlined in Sooter et al.
%
function [phic] = beta_lift(phic_reduc, beta)
phic_reduc = phic_reduc(:);
n = length(phic_reduc)+beta/2; %model order
A = zeros(beta/2);
b = zeros(beta/2, 1);
for i = 1:size(A, 1)
    A(i, :) = ((n-size(A, 2)+1):n).^(i-1);
    if i == 1
        b(i) = 1-sum(phic_reduc);
    else
        b(i) = -((1:length(phic_reduc)).^(i-1))*phic_reduc;
    end
end

phic = [phic_reduc ; A\b];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions
% checks whether the AR model with history kernel kern is explosive
%
%  i.e. Whether any roots of the AR polynomial are inside the unit circle
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization

function QC_plot(kern, kernc, db, b, deltaT)
figure; dim = length(kern)-b/2;
for i = 1:dim
    subplot(1, dim, i);
    hold on;
    xx = linspace(kernc(i)-0.2, kernc(i)+0.2, 100);
    klr = zeros(size(xx));
    expl = zeros(size(xx));
    for j = 1:length(xx)
        kernc_temp = kernc(1:dim); kernc_temp(i) = xx(j);
        kernc_temp = beta_lift(kernc_temp, b);
        klr(j) = log2(exp(1))/deltaT*get_klr_b(kern, kernc_temp(1:dim), b, 'trapz');
        expl(j) = explosive(kernc_temp);
    end
    plot(xx(expl==1), klr(expl==1), '.r', 'MarkerSize', 18);
    plot(xx(expl==0), klr(expl==0), '.c', 'MarkerSize', 18);

    xline(kernc(i), 'k', 'LineWidth', 2)
    yline(db, 'k', 'LineWidth', 2)
    set(gca, 'FontSize', 12);
    xlabel(['$(\phi_c)_' num2str(i) '$'], 'Interpreter', 'latex', 'FontSize', 17);
    ylabel('KL rate (bits/s)', 'FontSize', 14);
    set(gca, 'LineWidth', 1.5); box on;
    yl = ylim; ylim([0 yl(2)]);
end
end
