function [forcing] = act_forcing(c_R,c_G,c_E,t,...
                                 M,pts_eul,...
                                 R_sigma,G_sigma,E_sigma,...
                                 nu_RG,k_GE,...
                                 rate,...
                                 eta0,V_cn)
%Compute right hand side of ODE for (R*,G*,E*)
%   c_R,c_G,c_E are the nodal coefficients defining R*,G*,E* on mesh at
%    present time
%   t is the time at which we are evaluating
%   M = \int_vol phi_i*phi_j is the mesh mass matrix
%   n_pts is the number of nodes across the mesh
%   pts_eul is 3 x n_pts array giving Eulerian cone coordinates of nodes
%   R_sigma,G_sigma,E_sigma are the resp surface densities of all R,G,E
%   nu_RG is the activation rate of R for G
%   k_GE is the activation rate of G for E.
%   eta0 = .5*nu*epsilon_0
%   Vcone is the volume of the cone

%% R* Source Terms
F_1 = R_sigma;

% If flash is time dependent, source term contains extra term that could
%  not be lumped into K_gl
[~, flag_tdep] = flash(pts_eul,NaN,...
                       rate,true);

if flag_tdep
    F_1 = F_1 - c_R;
end

F_1 = F_1.*flash(pts_eul,t,rate);

F_1 = M*F_1;
f_R = eta0/(R_sigma*V_cn)*F_1;

%% G* Source Terms
F_2 = -nu_RG/G_sigma*((c_G+c_E).*c_R) + k_GE*(c_E.*c_G);

f_G = M*F_2;

%% E* Source Terms
F_3 = -k_GE*(c_E.*c_G);

f_E = k_GE*M*F_3;

%% Aggregate as Forcing Term
forcing = [f_R;...
           f_G;...
           f_E];
end

