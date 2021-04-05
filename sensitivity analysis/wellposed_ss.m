function [flag_proceed,flag_solvable] = wellposed_ss(datamat,pointer,...
                                                     M,Msl,Mfo)
%A priori check to see if parameters are consistent with flux balance
%   The nonlinear steady state equations lead to a system of constraints
%   between the parameters to see if they admit a solution for steady state
%   flux balance
%   flag_proceed is true or false beased on whether code may proceed with
%    this parameter set. This is decided by a combination of checking the
%    equations for their theoretical solvability and then checking if the
%    obtained steady state values fall in the ranges listed below.
%   flag_solvable indicates whether the given parameters admit a
%    theoretical steady state even if not physically admissible 
%   M, Msl, Mfo are optional args that the script will build if they are
%    not passed

%% Initialize and define acceptable physical ranges for quantities
flag_proceed = true;
flag_solvable = true;

% Accepted cG_d and Ca_d ranges
cG_d = [2 4]; % M
Ca_d = [.2 .29];% M (Korenbrot 2012 says .4)

% Accepted dark current ranges
j_d = [19.18 19.81]; % pA

%% Check acceptable ranges for parameters 
% Do this here to make sure that later operations are well-defined
if nnz(datamat<0) ~= 0
    % Some param isn't positive
    flag_proceed = false;
    flag_solvable = false;
    return
end

%if datamat(pointer(8)+3)>750
    % PDE_sigma too high
%    flag_proceed = false;
    %flag_solvable = false;
%    return
%end

%% Standardize the data
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
      
%% Check cases to see if solvable
% This steady state criterion is equivalent to the one in Shen
val1 = 1+(alpha_min/(Beta_dark*K_cG))^(m_cG);
val2 = (1-2*j_ex_sat/(j_cG_max*f_Ca))^(-1);
if val1 >= val2
    % Dark steady state does not exist 
    flag_proceed = false;
    flag_solvable = false;
    return
end

%% Check if cG_d and Ca_d values are permissible
% Compute steady state concentrations
[u_ss,v_ss]=steady_state(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,...
                                  flag_ch,...
                                  tol_stat,...
                                  B_Ca,F,...
                                  PDE_sigma,k_hyd,...
                                  u_tent,...
                                  v_tent,...
                                  j_cG_max,m_cG,K_cG,f_Ca,...
                                  j_ex_sat,K_ex,...
                                  alpha_max,alpha_min,m_cyc,K_cyc);


if (u_ss < cG_d(1))||(u_ss > cG_d(2))
    flag_proceed = false;
    return
elseif (v_ss < Ca_d(1))||(v_ss > Ca_d(2))
    flag_proceed = false;
    return
end

%% Check if dark current values are permissible
if nargin == 2
% Generate mesh
[pts,prisms,faces_sl,faces_fo,...
          n_pts,n_prism,n_fsl,n_ffo] = ...
          genmesh(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,sigma,...
                  n_sez,taglia,tol_R,tol_angle);
              
% Generate FEM matrices
M = mvol_ss(pts,prisms,...
         R_b,R_t,H);
Msl = mbd_ss(pts,faces_sl,...
          R_b,R_t,H);
Mfo = mbd_ss(pts,faces_fo,...
          R_b,R_t,H);
else
n_pts = size(M,1);
end

% Compute J_tot
Sigma_sl = sum(sum(Msl));
V_cn = sum(sum(M));
Sigma_fo = sum(sum(Mfo));

eta0 = .5*nu*epsilon_0;
Sigma_Ch = flag_ch(1)*Sigma_sl + ...
           flag_ch(2)*(nu/(1+nu))*V_cn/eta0 + ...
           flag_ch(3)*Sigma_fo;

% Expand u_ss and v_ss to full mesh
u_ss = repmat(u_ss,n_pts,1);
v_ss = repmat(v_ss,n_pts,1);
       
% cG-gated Current
J_cG = j_cG_max/Sigma_Ch*(u_ss.^m_cG)./(K_cG.^m_cG + u_ss.^m_cG);

% Exchanger Current
J_ex = j_ex_sat/Sigma_Ch*(v_ss./(K_ex+v_ss));

% Aggregate Current Density
J = J_cG + J_ex;

% Integrate over domain of channels
J_tot = flag_ch(1)*sum(Msl*J) + ...
           flag_ch(2)*(nu/(1+nu))/eta0*sum(M*J) + ...
           flag_ch(3)*sum(Msl*J);
       
% Normalize out the B_Ca*F factor from j_cg and j_ex
J_tot = B_Ca*F*J_tot;

% Convert units to pA
J_tot = (1e12)*J_tot;
J_tot = J_tot(1);

% Check if J_tot in right range
if (J_tot < j_d(1))||(J_tot > j_d(2))
    flag_proceed = false;
end

end

function [Mbd] = mbd_ss(pts,rects,...
                    R_b,R_t,H)
%Build (\int phi_i*phi_j dsigma) across dcone |_|rects

%% Needed parameters
n_pts = size(pts,2);
n_rects = size(rects,2);

%% Initialize lists used to build sparse M
% Number of nodal pairs we will cycle over while redundant (i,j) values
% will be accumulated by sparse
n_pairs = 16*n_rects;

I = zeros(n_pairs,1);
J = zeros(n_pairs,1);
V = zeros(n_pairs,1);

%% Loop over elements and build lists
for k=1:n_rects
    % Slice in needed values
    rect = rects(:,k);
    
    % Compute Mvol_loc
    Mbd_loc = mbd_loc(pts,rect,R_b,R_t,H);
    
    % Write these values in list format
    % Order of indices is compatible across I,J, V because list order is 
    % first down array column then onto next column
    
    %  Array value is i-index there
    ram_I = repmat([1;2;3;4],1,4);
    ram_I = ram_I(:);
    %  Now switch from local node indexing to mesh indexing
    ram_I = rect(ram_I);
    
    %  Array value is j-index there
    ram_J = repmat([1 2 3 4],4,1);
    ram_J = ram_J(:);
    %  Now switch from local node indexing to mesh indexing
    ram_J = rect(ram_J);
    
    ram_V = Mbd_loc(:);
    
    % Write these values into I,J,V
    I((k-1)*16+1:k*16) = ram_I;
    J((k-1)*16+1:k*16) = ram_J;
    V((k-1)*16+1:k*16) = ram_V;
end

%% Build Mbd
Mbd = sparse(I,J,V,n_pts,n_pts);

end

function [M] = mvol_ss(pts,prisms,...
                    R_b,R_t,H)
%Build (\int phi_i*phi_j dx) across cone |_|prisms

%% Needed parameters
n_pts = size(pts,2);
n_prisms = size(prisms,2);

%% Initialize lists used to build sparse M
% Number of nodal pairs we will cycle over while redundant (i,j) values
% will be accumulated by sparse
n_pairs = 36*n_prisms;

I = zeros(n_pairs,1);
J = zeros(n_pairs,1);
V = zeros(n_pairs,1);

%% Loop over elements and build lists
for k=1:n_prisms
    % Slice in needed values
    prism = prisms(:,k);
    
    % Compute Mvol_loc
    Mvol_loc = mvol_loc(pts,prism,R_b,R_t,H);
    
    % Write these values in list format
    % Order of indices is compatible across I,J, V because list order is 
    % first down array column then onto next column
    
    %  Array value is i-index there
    ram_I = repmat([1;2;3;4;5;6],1,6);
    ram_I = ram_I(:);
    %  Now switch from local node indexing to mesh indexing
    ram_I = prism(ram_I);
    
    %  Array value is j-index there
    ram_J = repmat([1 2 3 4 5 6],6,1);
    ram_J = ram_J(:);
    %  Now switch from local node indexing to mesh indexing
    ram_J = prism(ram_J);
    
    ram_V = Mvol_loc(:);
    
    % Write these values into I,J,V
    I((k-1)*36+1:k*36) = ram_I;
    J((k-1)*36+1:k*36) = ram_J;
    V((k-1)*36+1:k*36) = ram_V;
end

%% Build M
M = sparse(I,J,V,n_pts,n_pts);

end
