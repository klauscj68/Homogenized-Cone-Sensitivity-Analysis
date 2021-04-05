% Plot the total quantity of R*,G*,E* over the cone as measured by []_sigma

%% Load the data set
load('..\activation.mat');

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

%% Build M matrix for integrating
[pts,prisms,faces_sl,faces_fo,...
          n_pts,n_prism,n_fsl,n_ffo] = ...
          genmesh(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,sigma,...
                  n_sez,taglia,tol_R,tol_angle);
              
M = mvol(pts,prisms,...
                    R_b,R_t,H);
                
eta0 = .5*nu*epsilon_0;
                
%% Plot R* over time
Rst = sum(M*Rst_sig);

% Finish conversion from surface density to volumic
Rst = 1/eta0*Rst;
              
% Run plotting
figure
plot(taxis,Rst,'LineWidth',2);

title('Total Rst in Volume')

%% Plot G* Over Time
Gst = sum(M*Gst_sig);

% Finish conversion from surface density to volumic
Gst = 1/eta0*Gst;
              
% Run plotting
figure
plot(taxis,Gst,'LineWidth',2);

title('Total Gst in Volume')

%% Plot E* Over Time
Est = sum(M*Est_sig);

% Finish conversion from surface density to volumic
Est = 1/eta0*Est;
              
% Run plotting
figure
plot(taxis,Est,'LineWidth',2);

title('Total Est in Volume')

%% Reset path and vars
path(oldpath);
clear

