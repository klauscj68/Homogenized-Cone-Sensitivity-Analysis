function [k_hyd] = fk_hyd(nu,epsilon_0,Beta_dark,PDE_sigma)
%Compute k_hyd from its parameters

eta0 = .5*nu*epsilon_0;
k_hyd = eta0*Beta_dark/PDE_sigma;

end

