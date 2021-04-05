function [flag_cauchy] = cauchy_crit(Cauchy,...
                                     avg_Cauchy, norma_inf,...
                                     tol_fix)
%Test the samples to see if columns are sufficiently Cauchy
%   Cauchy is an array whose columns are the samples we wish to test for
%    Cauchy convergence.
%   avg_Cauchy is the number of columns we average together before testing
%    if Cauchy, to average out small numerical oscillation. In particular,
%    size(samples,2) must be divisible by cluster_sz
%   norma_inf = true means Linfty across node samples is used for
%    convergence. norma_inf = false means L1 across node samples is used.
%   tol_fix is how small the relative error between any two samples must be
%    to pass the convergence test 
%   flag_Cauchy is true if the samples pass the test and false if not

%% Initialize parameters
[n_pts,n_smp] = size(Cauchy);

% Number of averaged samples we directly apply Cauchy convergence too
tst_Cauchy = n_smp/avg_Cauchy;

%% Cluster average samples into their groups
% Average the candidate Cauchy samples in
%  blocks of size avg_Cauchy 
ram = zeros(n_pts,tst_Cauchy);
for j=1:tst_Cauchy
    ram(:,j) = mean(...
                    Cauchy(:,(j-1)*avg_Cauchy + 1:...
                              j*avg_Cauchy),...
                    2);
end

%% Restructure as successive differences
% Restructure ram(:,1:tst_Cauchy-1) as successive
%  differences
ram(:,1:tst_Cauchy-1) = ram(:,2:tst_Cauchy) - ...
                         ram(:,1:tst_Cauchy-1);
                     
%% Test for convergence
%  Set appropriate norm
if norma_inf == true
    % Error assessed by Linf norm
    p = Inf;
else
    % Error assessed by L1 norm
    p = 1;
end
                
%  Compute norms
err_Cauchy = zeros(1,tst_Cauchy);
for j=1:tst_Cauchy
    % Must apply norm column by column or else norm
    %  reverts to acting on full matrix, which is not
    %  wanted
    err_Cauchy(j) = norm(ram(:,j),p);
end
                
%  See if relative error is sufficiently small
%   By Tri-inequality, if full rel sum is underneath tol_R,
%   then so is norm between any of intermediate samples
r_error = sum(err_Cauchy(1:tst_Cauchy-1));

if (r_error <= tol_fix*err_Cauchy(:,tst_Cauchy))||...
   ( (err_Cauchy(:,tst_Cauchy) <= 9e-16)&&...
      r_error <= min(1e-13,tol_fix) )
    % Passed the test
    flag_cauchy = true;
else
    % Failed test
    flag_cauchy = false;
end

end

