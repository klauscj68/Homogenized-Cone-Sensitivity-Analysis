function [F,G,H,K,...
          mu,sigma,arms,...
          nqty] = mc_plot(qty,n_smp,threshold)
%Generate an empirical plot of sample avg of qty evidencing MC convergence
%   qty is a 1 x n_trials vector recording the value of integrand in each
%    trial
%   n_smp is the number of points we wish to have in our plots.  Roughly we
%   partition qty into consecutive n_binsz = n_trials/n_smp sized bins and 
%   average the values there to become a single sample for that bin
%   F,G are the object handles to the mean and normal distribution plots
%   This convergence check is based on the argument in Ch 3 of Monte Carlo
%    Statistical Methods p. 84
%   mu and sigma are the sample mean and standard deviation for ALL the
%    samples, not after they've been clustered into bins of size 100. 
%   arms stands for average root mean square error for the qq plot between
%    sample values and the line y = x
%   nqty is the number of elements in qty taken as convergent
%% Declare whether generating figures
flag_plot = false;

%% Optionally exclude identified nonconvergent runs (default is dont)
thresh = 0; % units of percent
% Override if thresh value was passed
if nargin == 3
    thresh = threshold;
end
rawqty = qty;
[stqty,idx] = sort(rawqty);
nqty = length(rawqty);

if mod(nqty,2) == 0
    middle = .5*( stqty(nqty/2)+stqty(nqty/2+1) );
else
    middle = stqty( (nqty-1)/2 + 1 );
end

nexcl = nqty*thresh*1e-2;
nexcl = floor(nexcl);

% Truncate the data
%  Find the values furthest from the median
outlier = stqty - middle;
[~,jdx] = sort(abs(outlier));

%  Reduce to chosen set
qty( idx(jdx(end-nexcl+1:end)) ) = [];
nqty = length(qty);
%% Output the TOTAL sample mean and standard deviation
mu = mean(qty);
sigma = sqrt(...
             mean( ...
                  (qty-mu).^2 ...
                 ) ...
            );
%% Generate histogram of samples
if flag_plot
    K = figure;
    histogram(qty,'Normalization','probability');
    title(['Data: ' num2str(length(qty)) ' samples']);
else
    K = NaN;
end
%% Initialize
n_trials = max(size(qty));
n_binsz = floor(n_trials/n_smp);

if n_binsz == 0
    % Number of desired smp exceeded number of trials
    n_smp = n_trials;
    n_binsz = 1;
end

% Truncate trials to be evenly divisible by n_binsz
n_trials = n_smp*n_binsz;
try
    qty = qty(1:n_trials);
catch
    F = NaN; G = NaN; H = NaN; K = NaN;
    mu = NaN; sigma = NaN; arms = NaN;
    return
end
% Reshape qty into n_binsz by n_smp
%  Each column is to be a sample point in plots
qty = reshape(qty,n_binsz,n_smp);

%% Mean plot
smp_tot   = n_binsz*(1:n_smp);
smp_means = cumsum(sum(qty,1));
smp_means = smp_means./smp_tot;

if flag_plot
    F = figure;
    hold on
    plot(smp_tot,smp_means,'b','LineWidth',3);
    plot(smp_tot,smp_means(n_smp)*ones(1,n_smp),'r','LineWidth',3)

    F.CurrentAxes.FontSize = 16;
    F.CurrentAxes.TitleFontSizeMultiplier = 1.25;
    xlabel('Number of samples');
    ylabel('Sample of Means');
    title('MC Mean Convergence')
else
    F = NaN;
end

%% Normal Distr Plot
% Compute the sample means which are ~ N(mu,sigma/\sqrt{N})
smp_means = mean(qty,1);

% The approximating normal distribution parameters
avg = mean(smp_means);
vari = mean( (smp_means - avg).^2 );

% Normalize to standard normal
Nsmp_means = (smp_means - avg)/sqrt(vari);

if flag_plot
    G = figure;
    hold on
    histogram(Nsmp_means,'Normalization','pdf');

    NPDF = @(x) (1./sqrt(2*pi)).*exp(-(x.^2)./2);
    fplot(NPDF, G.CurrentAxes.Children(1).BinLimits,...
            'r','LineWidth',3);

    G.CurrentAxes.FontSize = 16;
    G.CurrentAxes.TitleFontSizeMultiplier = 1.25;
    xlabel('Value Normalized Sample Mean');
    ylabel('PDF');
    title('MC Normalized CLT Convergence')
else
    G = NaN;
end

%% Q-Q Plot
% Compare the sample means with a N(0,1) normal distribution to evidence
%  convergence
Nsmp_means_ord = unique(Nsmp_means);

for i=length(Nsmp_means_ord):-1:1
    CDF(:,i) = [nnz(Nsmp_means <= Nsmp_means_ord(i))/length(Nsmp_means);...
                normcdf(Nsmp_means_ord(i))];
end
CDF = CDF';

% Setup the xy line
xyline = linspace(0,1,1e3);

% Compute the arms erros
arms = sqrt(mean( (CDF(:,2)-CDF(:,1)).^2 ));

% Make the plot
if flag_plot
    H = figure;
    hold on

    plot(CDF(:,1),CDF(:,2),'bx')
    plot(xyline,xyline,'r','LineWidth',3);

    H.CurrentAxes.FontSize = 16;
    H.CurrentAxes.TitleFontSizeMultiplier = 1.25;
    xlabel(['Empirical CDF: arms =' num2str(arms,4)]);
    ylabel('Standard Normal CDF');
    title('Q-Q Plot of Sample Means')
else
    H = NaN;
end
end

