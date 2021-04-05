function [F] = plot_sliver(u,taxis,t,...
                              var_name)
%Plot full lateral boundary of cone and colormap the sliver section
%   u is the n_pts x n_tpts solution array
%   t is the time point at which we are sampling
%   varname is the string with the variable being plotted

%% Set path and get data
oldpath = addpath('..\common\');
addpath('..\ode solver\');

[datamat,pointer,...
          fnames] = data_set();

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
      
%% Regenerate mesh
[pts,prism,faces_sl,faces_fo,...
          n_pts,n_prism,n_fsl,n_ffo] = ...
          genmesh(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,sigma,...
                  n_sez,taglia,tol_R,tol_angle);
              
%% Map points to eulerian cone
pts_eul = map_to_eul(pts,R_b,R_t,H);

%% Sample u nodals at time point
u = sol_interp(u,taxis,t);

%% Plot the colormap
F = figure; hold on
for i=1:n_fsl
    % Geometric locations of face nodes in this rectangle
    xvert = pts_eul(1,faces_sl(:,i)); xvert = xvert(:);
    yvert = pts_eul(2,faces_sl(:,i)); yvert = yvert(:);
    zvert = pts_eul(3,faces_sl(:,i)); zvert = zvert(:);
    
    % Fill 3 needs these things listed in oriented order, not with twist
    %  like was the convention when I wrote the faces lists
    P = eye(4);
    P = [P(:,1) P(:,2) P(:,4) P(:,3)];
    xvert = P*xvert;
    yvert = P*yvert;
    zvert = P*zvert;
    
    % Nodal values of u at these locations
    cval = u(faces_sl(:,i)); cval = cval(:);
    
    % Plot the patch
    fill3(xvert,yvert,zvert,cval,'LineStyle','none');
end

% Optional segment to flip colors upside down if too much of sliver is in
%  yellow
c = parula;
c = flipud(c);
colormap(c);

% Return to normal plotting
colorbar;
F.CurrentAxes.FontSize = 20;
F.CurrentAxes.FontName = 'TimesNewRoman';
t_round = floor(t*1000); t_round = t_round/1000;
title([var_name ' at time ' num2str(t_round)])
view(3);

%% Return path to normal
path(oldpath);
end

