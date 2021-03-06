% Simulate the cascade for light activation of cG and Ca
%% Grab the data
        
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
          data(datamat,pointer);
      
%% Generate the mesh
[pts,prisms,faces_sl,faces_fo,...
          n_pts,n_prism,n_fsl,n_ffo] = ...
          genmesh(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,sigma,...
                  n_sez,taglia,tol_R,tol_angle);
              
pts_eul = map_to_eul(pts,...
                     R_b,R_t,H);
      
%% Generate Volumetric FEM matrices
M = mvol(pts,prisms,...
         R_b,R_t,H);
K = kvol(pts,prisms,...
         R_b,R_t,H);
Msl = mbd(pts,faces_sl,...
          R_b,R_t,H);
Ksl = kbd(pts,faces_sl,...
          R_b,R_t,H);
Mfo = mbd(pts,faces_fo,...
          R_b,R_t,H);

%% Integrate the FEM System
% Compute steady state values
[cG0,Ca0]=steady_state(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,...
                       flag_ch,...
                       tol_stat,...
                       B_Ca,F,...
                       PDE_sigma,k_hyd,...
                       u_tent,...
                       v_tent,...
                       j_cG_max,m_cG,K_cG,f_Ca,...
                       j_ex_sat,K_ex,...
                       alpha_max,alpha_min,m_cyc,K_cyc);

S0_gl = [repmat(cG0,n_pts,1);...
         repmat(Ca0,n_pts,1)];


% Initialize remaining parameters for activation
load('Evol.mat');
E_st = Evol;
tpts = taxis;
flag_case = 'Msg';
     
% Run integration
[S_gl,taxis] = ode_integrate(M,K,Msl,Ksl,Mfo,n_pts,...
                                      pts_eul,...
                                      S0_gl,...
                                      E_st,tpts,...
                                      flag_case,...
                                      flag_preassembled,...
                                      R_b,R_t,H,cosgamma0,theta_in,theta_fin,epsilon_0,nu,sigma,...
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
                                      flag_restart,rst_index);
                                  
%% Break solution up into cG,Ca
cG = S_gl(1:n_pts,:);
Ca = S_gl(n_pts+1:2*n_pts,:);

%% Compute J_tot and J_drop
Sigma_sl = sum(sum(Msl));
V_cn = sum(sum(M));
Sigma_fo = sum(sum(Mfo));

eta0 = .5*nu*epsilon_0;
Sigma_Ch = flag_ch(1)*Sigma_sl + ...
           flag_ch(2)*(nu/(1+nu))*V_cn/eta0 + ...
           flag_ch(3)*Sigma_fo;
       
% cG-gated Current
J_cG = j_cG_max/Sigma_Ch*(cG.^m_cG)./(K_cG.^m_cG + cG.^m_cG);

% Exchanger Current
J_ex = j_ex_sat/Sigma_Ch*(Ca./(K_ex+Ca));

% Integrate Current Density
J = J_cG + J_ex;

% Integrate over domain of channels
J_tot = flag_ch(1)*sum(Msl*J) + ...
           flag_ch(2)*(nu/(1+nu))/eta0*sum(M*J) + ...
           flag_ch(3)*sum(Msl*J);

J_drop = 100*(J_tot(1) - J_tot)/J_tot(1);

% Normalize out the B_Ca*F factor from j_cg and j_ex
J_tot = B_Ca*F*J_tot;


save('messenger','cG','Ca','taxis',...
     'J_tot','J_drop','-v7.3')
