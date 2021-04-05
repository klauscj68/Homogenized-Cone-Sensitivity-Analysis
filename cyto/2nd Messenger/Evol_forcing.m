function [forcing] = Evol_forcing(c_E,...
                                  G_star,...
                                  M,...
                                  E_vol,...
                                  k_GE)
%Compute right hand side of ODE for E_vol
%   c_E are the nodal coefficients of E in column form at present time
%   G_star is the nodal coefficients of G_star at present time
%   M is the \int phi_i*phi_j matrix over the volume mesh
%   n_pts is the number of nodes in the mesh
%   E_vol is the total volumic concentration of all E in the cone
%   k_GE is the activation rate of E for G

%% Activation of E* by G*
forcing = k_GE*(E_vol - c_E).*G_star;
forcing = M*forcing;

end

