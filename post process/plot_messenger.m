% Plot total quantity of cG and Ca2+ over cone and their current (drops)
%% Load the data set
load('..\messenger.mat');

%% Add needed files to path
oldpath = addpath('..\common\');
addpath('..\cyto\');

%% Grab data parameters
[datamat,pointer] = data_set();
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
          alpha_max,alpha_min,m_cyc,K_cyc] = ...
          data(datamat,pointer);

%% Build M matrices for integrating
[pts,prisms,faces_sl,faces_fo,...
          n_pts,n_prism,n_fsl,n_ffo] = ...
          genmesh(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,sigma,...
                  n_sez,taglia,tol_R,tol_angle);
              
M = mvol(pts,prisms,...
                    R_b,R_t,H);
                
M_sl = mbd(pts,faces_sl,...
                        R_b,R_t,H);
                    
M_fo = mbd(pts,faces_fo,...
                        R_b,R_t,H);
                    
%% Volume Normalization
V_cn     = sum(sum(M));
                
%% Plot average [cG] over time
cG_avg = sum(M*cG)/V_cn;
              
% Run plotting
figure
plot(taxis,cG_avg,'LineWidth',2);

title('Avg [cG] in Volume')

%% Plot [Ca] over time
Ca_avg = sum(M*Ca)/V_cn;
              
% Run plotting
figure
plot(taxis,Ca_avg,'LineWidth',2);

title('Avg [Ca] in Volume')

%% Plot Current over Time
figure
% Normalize J_tot to pA and x -1
J_tot = -J_tot;
J_tot = J_tot*10^(12);
plot(taxis,J_tot - J_tot(1),'LineWidth',2);
ylabel(['-\Delta' 'j (pA)']);
title('Total Current over Channels')

%% Plot Drop over Time
figure

plot(taxis,J_drop,'LineWidth',2);

title('Current Drop over Time')
ylabel('% of J_d_a_r_k')

%% Reset path and vars
path(oldpath);
clear