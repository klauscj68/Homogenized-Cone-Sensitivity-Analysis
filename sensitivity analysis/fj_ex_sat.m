function [j_ex_sat_normalized] = fj_ex_sat(j_ex_sat_normalized,B_Ca,F,...
                                           delta_r)
%Compute j_ex_sat from its parameters
%   See discussion in fj_cG_max

B_Ca_old = B_Ca/(1+delta_r);
j_ex_sat = j_ex_sat_normalized*B_Ca_old*F;

j_ex_sat_normalized = j_ex_sat/(B_Ca*F);
end