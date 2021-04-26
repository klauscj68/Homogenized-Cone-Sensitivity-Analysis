function alpha_min = falpha_min(amin,alpha_max,Beta_dark,K_cG,j_ex_sat,j_cG_max,f_Ca,m_cG)
%Compute alpha_min using amin, alpha_max and other params
%   A valid choice of alpha_min must satisfy the steady-state criterion and
%    be beneath alpha_max. Here alpha_min is computed as a function of
%    these choices

alpha_min = amin*min(alpha_max,...
                     Alimit(Beta_dark,K_cG,j_ex_sat,j_cG_max,f_Ca,m_cG));

end

function Amax = Alimit(Beta_dark,K_cG,j_ex_sat,j_cG_max,f_Ca,m_cG)
% Return max value of alpha_min permitted for well-posed steady state 
%  Note: formula safeguards against errors in j_ex_sat,j_cG_max
%  normalization by 1/B_Ca*F since its the ratio of currents that is used
    Amax = Beta_dark*K_cG*...
               ((1-2*j_ex_sat/(j_cG_max*f_Ca))^(-1) - 1)^(1/m_cG);
end
