
% Add path to actual toolbox
addpath '../src/'



%%%%%% Sam's Demo (from demo.mlx) %%%%%%%%%%%%%%%%%%%%
if (0)
close all;
kern = [0.8 0.2 -0.1]; sigma = 1; % AR model parameters
T = 1e4; % number of time steps
x = runAR(kern, sigma, T); % simulate AR model
end

%%%%% Generating synthetic data %%%%%%%%%%
%close all;
n = 12; p = 0.7;
if n >6
    kern = (-p).^([0:n-1]);   sigma = 1; % AR model parameters
else
    kern = (p).^([0:n-1]);   sigma = 1; 
end
kern = 0.95*kern/sum(kern);

% However I make the kernel, it's possible I come up with 
%  an explosive one. Scale down kernel until I get something ok
minroot = min(abs(roots([fliplr(kern(:)') -1])));
while minroot < 1
    kern = kern * minroot * 0.99;
    minroot = min(abs(roots([fliplr(kern(:)') -1])));
    disp('Adjusting kernel for stability')

end

T = 1e4; % number of time steps
x = runAR(kern, sigma, T); % simulate AR model


% d2 settings (minimal)
order = length(kern); % AR model order
deltaT = 1; % time between consecutive steps in time series (in s)
b = 6; % calculate d2 (rather than d4, d6, ...)
with_err_bars = false; % don't calculate error bars
with_QC = false; % don't make quality control plots
with_parallel = false; % don't use parallel computing
fit_method = 'YuleWalker'; % AR fit method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep my original kernel and sigma
%
% i.e. return "kernf" and "sigmaf" as the
%  fitted values

% Note first step is to normalize x to have unit variance

[db, sddb, kernf, sigmaf, H, kernc, exit_status] = calc_db(x, order, deltaT, ...
    b, with_err_bars, with_QC, with_parallel, fit_method);

disp(['original kernel: [' num2str(kern(:)') ']'])

disp(['best-fit kernel: [' num2str(kernf(:)') ']'])

disp(['closest to best-fit kernel on Mb: [' num2str(kernc(:)') ']'])


%% Compare original kernal used to generate synthetic data, 
%% fitted kernel, and the kernel on the critical manifold
figure;
subplot(2,1,1);
% True kernel
plot(kern,'*'); hold on;
plot(kernf, '*');  plot(kernc, '*');
set(gca,'FontSize',16); legend('True kernel', 'Best-fit','Best fit, critical')
subplot(2,1,2);
% AR data
plot(x); set(gca,'FontSize',16);


%%%%%%%%%%%% AUXILIARY %%%%%%%%%%%%%%%%
function x = runAR(kern, sigma, T)
kern = kern(:);
x = zeros(1, T);
for t = 1:T
    hist = fliplr(x(max(1, t-length(kern)):(t-1)));
    x(t) = hist*kern(1:length(hist)) + sigma*randn;
end
end

function x = get_spike_counts(st, t1, t2, deltaT)
    bins = t1:deltaT:t2;
    x = histc(st, bins);
end