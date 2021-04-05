function [M_gl,K_gl] = Evol_assembly(M,K,...
                                     D_E_st,k_E)
%Build the matrices used in solving FEM system for volumic E
%   M = \int phi_i*phi_j and K = \int D_x,y phi_i* D_x,y phi_j are the
%       geometric matrices
%   The rest are biochemical parameters for E

%% Build M_gl
M_gl = M;

%% Build K_gl
K_gl = D_E_st*K + k_E*M;

end

