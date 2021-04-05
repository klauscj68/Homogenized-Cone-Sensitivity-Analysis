function [value,functional] = fitness(index)
%Compute fitness of given data values for matching desired features
%   It's assumed the data is saved in trial_i.mat where i is given by
%   index. The function then computes fitness by adding together by what
%   relative percentage of the desired interval width is the parameter
%   values outside the interval
%   Desired
%   rmax = 21\pm .9 pA
%   I_1/2 = 940 
%   t_peak \in [33,59]
%   [cG}_d = [2-4] \muM
%   [Ca^{2+}]_d = [200-670] nM
%   rmax = [.02,.05]*max current

%% Initialize
functional = zeros(7,1);

%% Load the data
load(['trial_' num2str(index) '.mat']);

[R_b,R_t,H,cosgamma0,theta_in,theta_fin,epsilon_0,nu,sigma,...
          flag_ch,...
          n_sez,taglia,tol_R,tol_angle,...
          method_cyto,theta,alpha,tol_fix,norma_inf,t_fin,...
          n_step_t,downsample,...
          tol_stat,...
          N_Av,F,...
          B_cG,B_Ca,...
          nu_RG,D_R_st,k_R,R_sigma,...
          k_GE,D_G_st,G_sigma,...
          D_E_st,k_E,PDE_sigma,Beta_dark,K_m,k_cat,...
          E_sigma,E_vol,k_hyd,kcat_DIV_Km,k_st,...
          u_tent,kk_u,...
          v_tent,kk_v,...
          j_cG_max,m_cG,K_cG,f_Ca,...
          j_ex_sat,K_ex,...
          alpha_max,alpha_min,m_cyc,K_cyc,...
          n_Rst0,rate,...
          flag_restart,rst_index] = ...
          data(smp,pointer);
      
%% Experimental trace l2 error
% Reason for these value choices is the output of flashimg
%{
%# Fitting for Fain Mouse
%  Interval is range on coordinates of flashimg
t_interval = [0.0177 0.35];
%  Experimental trace placed flash at 5ms
tshift = .005; 

% Broke the image import up into three subintervals
%  bc fitting routine was more robust this way
p1 = @(x) (-2.509e5)*x.^3+(1.563e4)*x.^2+206.5*x-5.33;
p2 = @(x) (141.3)*x.^2-194.9*x+21.74;
p3 = @(x) -4253*x.^3+2262*x.^2-400.2*x+23.31;

flashfit = @(x) p1(x).*(x<=0.046699) + ...
                p2(x).*(x> 0.046699).*(x<=0.1185) +...
                p3(x).*(x>0.1185).*(x<=0.19);
            
% Compute l2 error of this trial with polynomial
%  Find times in coordinates of model that are in domain of fit
flag_t = logical((taxis >= (t_interval(1)-tshift)).*...
                 (taxis <= (t_interval(2))-tshift));

% Note experimental trace placed a flash at 5 ms
tpts = taxis(flag_t);
%  To evaluate fit go back to domain of fit
tpts = tpts + tshift;
%}

%# Fitting for Korenbrot striped bass
%  Interval is range on coordinates of flashimg
t_interval = [0.001 0.25];

% Broke the image import up into three subintervals
%  bc fitting routine was more robust this way
p1 = @(x) 1.499e+05*x.^4+(-3.868e+04)*x.^3+1192*x.^2+255*x-0.8726;
p2 = @(x) 2.834e+04*x.^4-8055*x.^3-766*x.^2+241*x+1.783;
p3 = @(x) (1.088283552608914e+07)*x.^5+(-1.271803803813630e+07)*x.^4+(5.896281563293114e+06)*x.^3 ...
           +(-1.354415750606520e+06)*x.^2+(1.539943408224258e+05)*x+(-6.925865471730017e+03);

flashfit = @(x) p1(x).*(x<=0.1033) + ...
                p2(x).*(x>0.1033).*(x<=0.1987) +...
                p3(x).*(x>0.1987).*(x<=0.2504);
            
% Compute l2 error of this trial with polynomial
%  Find times in coordinates of model that are in domain of fit
flag_t = logical((taxis >= (t_interval(1))).*...
                 (taxis <= (t_interval(2))));

tpts = taxis(flag_t);

%# Resuming standard calculuation

simpts = J_tot(1)-J_tot; 
simpts = (1e12)*simpts;
simpts = simpts(flag_t);

flashpts = flashfit(tpts);

functional(1) = sqrt(...
                     sum(...
                         (flashpts - simpts).^2 ...
                        )...
                    );

%% Dark current value
interval = [21-.9,21+.9];
functional(2) = re(J_tot(1)*10^12,interval);

%% Half-maximal value
interval = [45 55];
functional(3) = re(max(J_drop),interval);

%% Time to peak drop value
[~,t_peak] = max(J_drop);
t_peak = taxis(t_peak);

interval = [33e-3 59e-3];
functional(4) = re(t_peak,interval);

%% cG dark value
interval = [2 4];
functional(5) = re(cG0,interval);

%% Ca dark value
interval = [200e-3 670e-3];
functional(6) = re(Ca0,interval);

%% Dark current percentage of total
% Magnitude of current and scale out the B_Ca*F factor dataset had included
% in
J_mag = j_cG_max + j_ex_sat;
J_mag = J_mag*B_Ca*F;

interval = [.02 .05];
functional(7) = re(J_tot(1)/J_mag,interval);

%% Compute desired error functional
value = functional(1);

end

function [rel_err] = re(point,interval)
%Compute the magnitude relative error point is in or outside interval
%   point is a single number
%   interval is a 1 x 2 giving left and right end points of intervals

% Check interval ordering
if interval(2) < interval(1)
    interval = [interval(2) interval(1)];
end

% Check where point falls
check = (point <= interval);
check = sum(check);

% Get interval width
delta = interval(2) - interval(1);

% Compute rel_err
switch check
    case 2
        % Point is beneath full interval
        rel_err = 100*(interval(1) - point)/delta;
        
    case 1
        % Point is in the full interval
        rel_err = 0;
        
    case 0
        % Point is above the full interval
        rel_err = 100*(point - interval(2))/delta;
        
end
end