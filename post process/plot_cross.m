function [u_avg] = plot_cross(u,taxis,t,...
                              z,...
                              var_name)
%Plot colormap of 3D solution at local or averaged cross section
%   u is the n_pts x n_tpts solution array
%   t is the time point at which we are sampling
%   pts, prism, taxis are the mesh for space and time
%   z should either be empty, 1 x 1, or 1 x 2.  If empty the code averages
%    across the entire cone.  If 1 x 1, the code shows the nearest z-height
%    cross section, and if 1 x 2 it averages over the cross sections of the
%    mesh in that range.
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


% Number of points in cross section
n_cpts = n_pts/n_sez;

% The product factor geometries
%  Return the unique z's and the first corresponding column of pts
[zaxis,jpts] = unique(pts(3,:));
cpts = pts(1:2,1:n_cpts);

% Extract triangles as first three rows of prism
n_tri = n_prism/(n_sez - 1);
tri = prism(1:3,1:n_tri);

%% Check that code expectations of data are true
assert(size(zaxis,2) == n_sez,'axis didnt have correct number sections');

check = jpts(2:n_sez) - jpts(1:n_sez-1);
check = unique(check);
assert( (size(check,1) == 1) && (size(check,2) == 1),...
        'space btwn sections wasnt even');
assert(check == n_cpts,'space btwn sections wasnt n_cpts');


%% Convert all cases to the interval case
if isempty(z)
    z = [0 max(pts(3,:))];
elseif (size(z,1) == 1) && (size(z,2) == 1)
    z = repmat(z,1,2);
end

%% Identify spanned column range of pts for given z range
i1 = find(z(1) >= zaxis,1,'last');
i1 = jpts(i1);
        
i2 = find(zaxis >= z(2),1);
i2 = jpts(i2);
% Since jpts gave first column with this z - value we need to pad columns
%  with rest of this section
i2 = i2 + (n_cpts - 1);

% Compute number of present sections
n_actsez = (i2 - (i1-1))/n_cpts;
assert(n_actsez == floor(n_actsez),'Did not find complete sections');

%% Sample u nodals at time point
u = sol_interp(u,taxis,t);

%% Accumulate and average nodals in this section
u_avg = u(i1:i2);
u_avg = reshape(u_avg,n_cpts,n_actsez);
u_avg = mean(u_avg,2);

%% Plot the colormap
F = figure; hold on;
for i=1:n_tri
    x_tri = cpts(1,tri(:,i)); x_tri = x_tri(:);
    y_tri = cpts(2,tri(:,i)); y_tri = y_tri(:);
    v_tri = u_avg(tri(:,i)); v_tri = v_tri(:);
    
    fill(x_tri,y_tri,v_tri,'LineStyle','none')
end
axis tight
colorbar;
F.CurrentAxes.FontSize = 20;
F.CurrentAxes.FontName = 'TimesNewRoman';
t_round = floor(t*100); t_round = t_round/100;
title([var_name ' at time ' num2str(t_round)])

%% Return path to normal
path(oldpath);

end

