function [E_vol] = fE_vol(PDE_sigma,nu,epsilon_0)
% Compute E_vol from its parameters

E_vol = 2*PDE_sigma/...
        (.5*nu*epsilon_0);

end

