function [M_gl,K_gl] = act_assembly(M,K,n_pts,...
                                    pts_eul,...
                                    D_R_st,D_G_st,D_E_st,...
                                    k_R,k_E,...
                                    R_sigma,G_sigma,E_sigma,...
                                    nu_RG,k_GE,...
                                    rate,...
                                    eta0,V_cn)
%Build the matrices used in solving FEM system for surface (R*,G*,E*)
%   Given the geometric M = \int_vol phi_i*phi_j and 
%   K = \int_vol D_xy phi_i*D_xy phi_j matrices assembled over the cone
%   volume, build the global matrices M_gl, K_gl for our FEM discretized
%   (R*,G*,E*) system.
%   n_pts is the number of nodes over the mesh
%   pts_eul are the nodes over Eulerian cone mesh
%   D's are the diffusion coefficients
%   k_R is the decay constant for activated R
%   k_E is the decay constant for the G-protein, Effector complex
%   eta0 = .5*nu*epsilon_0
%   Vcone is the volume of the cone

%% Build M_gl
% This matrix has diagonal block form
% [[M 0 0];...
%  [0 M 0];...
%  [0 0 M]]
[I_M,J_M,V_M] = find(M);

% Extract indices and values for sparse M_gl form
I_gl = [I_M;...
        n_pts + I_M;...
        2*n_pts + I_M];
J_gl = [J_M;...
        n_pts + J_M;...
        2*n_pts + J_M];
V_gl = repmat(V_M,3,1);

% Build Matrix
M_gl = sparse(I_gl,J_gl,V_gl);

%% Build K_gl
% Check whether the flash is ind of time
[~, flag_tdep] = flash(pts_eul,NaN,...
                       rate,true);

% If flash dependent on time our matrix is has diagonal block form
% [[D_R_st*K + k_R*M   0                                         0];
%  [-nu_RG*G_sigma/G_sigma*M  D_G_st*K + k_GE*E_sigma*M                 0];
%  [0                         -k_GE*E_sigma*M            D_E_st*K + k_E*M]]
        
[I_K,J_K,V_K] = find(K);

% Extract indices and values for sparse K_gl form, appealing to accumarray 
% feature of sparse

I_gl = [I_K;... (1,1)
        I_M;...
        n_pts + I_M;... (2,1)
        n_pts + I_K;... (2,2)
        n_pts + I_M;...
        2*n_pts + I_M;... (3,2)
        2*n_pts + I_K;... (3,3)
        2*n_pts + I_M];
    
J_gl = [J_K;... (1,1)
        J_M;...
        J_M;... (2,1)
        n_pts + J_K;... (2,2)
        n_pts + J_M;... 
        n_pts + J_M;... (3,2)
        2*n_pts + J_K;... (3,3)
        2*n_pts + J_M];

V_gl = [D_R_st*V_K;... (1,1)
        k_R*V_M;...
        -nu_RG*V_M;... (2,1)
        D_G_st*V_K;... (2,2)
        k_GE*E_sigma*V_M;...
        -k_GE*E_sigma*V_M;... (3,2)
        D_E_st*V_K;... (3,3)
        k_E*V_M];

        
if flag_tdep == false
% This matrix contains in its (1,1) block an additional summand
% eta0/(Vcone*R_sigma)*{Mflash.*(-)} 
%   This is of course equivalent to multiplying by a diagonal matrix

% Sample the time independent flash
rate = flash(pts_eul,NaN,...
             rate);
         
 % Convert to M[flash.*(-)] operator form
rate = spdiags(rate,0,size(M,1),size(M,2));
rate = eta0/(V_cn*R_sigma)*rate;
[I_flash,J_flash,V_flash] = find(M*rate);

I_gl = [I_gl;...
        I_flash];

J_gl = [J_gl;...
        J_flash];
    
V_gl = [V_gl;...
        V_flash];
    
end
        
% Build Matrix
K_gl = sparse(I_gl,J_gl,V_gl);

end

