function [R_b,R_t,H,cosgamma0,theta_in,theta_fin,epsilon_0,nu,sigma,...
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
          data(datamat,pointer)
%Extract array of datamat into its variable names for use in subroutines

%% Geometric parameters
R_b = datamat(1);
R_t = datamat(2);
H   = datamat(3);
cosgamma0 = datamat(4);
theta_in = datamat(5);
theta_fin = datamat(6);
epsilon_0 = datamat(7);
nu = datamat(8);
sigma = datamat(9);
flag_ch = logical(...
          datamat(10:12)...
          );
      
%% Space discretization
pos = pointer(1);
n_sez = datamat(pos + 1);
taglia = datamat(pos + 2);
tol_R = datamat(pos + 3);
tol_angle = datamat(pos + 4);

%% Time integration
pos     = pointer(2);
method_cyto = datamat(pos + 1: pos + 3);
theta   = datamat(pos + 4);
alpha   = datamat(pos + 5);
tol_fix = datamat(pos + 6);
norma_inf = datamat(pos + 7);
t_fin = datamat(pos + 8);
n_step_t = datamat(pos + 9);
downsample = datamat(pos+10);

%% Steady state solver
pos = pointer(3);
tol_stat = datamat(pos + 1);

%% Physical constants
pos = pointer(4);
N_Av = datamat(pos+1);
F = datamat(pos+2);

%% Cytoplasmic buffers
pos = pointer(5);
B_cG = datamat(pos+1);
B_Ca = datamat(pos+2);

%% R_st biochemistry
pos = pointer(6);
nu_RG = datamat(pos+1);
D_R_st = datamat(pos+2);
k_R = datamat(pos+3);
R_sigma= datamat(pos+4);

%% G_st biochemistry
pos = pointer(7);
k_GE = datamat(pos+1);
D_G_st = datamat(pos+2);
G_sigma = datamat(pos+3);

%% PDE_st biochemistry (free parameters)
pos = pointer(8);
D_E_st = datamat(pos+1);
k_E = datamat(pos+2);
PDE_sigma = datamat(pos+3);
Beta_dark = datamat(pos+4);
K_m = datamat(pos+5);
k_cat = datamat(pos+6);

%% PDE_st biochemistry (derived parameters)
pos = pointer(9);
E_sigma = datamat(pos+1);
E_vol = datamat(pos+2);
k_hyd = datamat(pos+3);
kcat_DIV_Km = datamat(pos+4);
k_st = datamat(pos+5);

%% cG biochemistry
pos = pointer(10);
u_tent = datamat(pos+1);
kk_u = datamat(pos+2);

%% Ca biochemistry
pos = pointer(11);
v_tent = datamat(pos+1);
kk_v = datamat(pos+2);

%% cG-gated current biochemistry
pos = pointer(12);
j_cG_max = datamat(pos+1);
m_cG = datamat(pos+2);
K_cG = datamat(pos+3);
f_Ca = datamat(pos+4);

%% Exchanger current biochemistry
pos = pointer(13);
j_ex_sat = datamat(pos+1);
K_ex = datamat(pos+2);

%% Cyclase biochemistry
pos = pointer(14);
alpha_max = datamat(pos+1);
alpha_min = datamat(pos+2);
m_cyc= datamat(pos+3);
K_cyc = datamat(pos+4);

%% Illumination parameters
pos = pointer(15);
n_Rst0 = datamat(pos+1);
rate = datamat(pos+2);

%% ODE restart parameters
pos = pointer(16);
flag_restart = logical(datamat(pos+1));
rst_index = datamat(pos+2);

end

