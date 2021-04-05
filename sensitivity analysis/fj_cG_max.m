function [j_cG_max_normalized] = fj_cG_max(j_cG_max_normalized,B_Ca,F,...
                                           delta_r)
%Compute j_cG_max from its parameters
%   j_cG_max_normalized = j_cG_max/(B_Ca*F)
%   Because of this subtlety we must first return j_cG_max to its
%   unnormalized form and then scale
%   delta_r is the amount which B_Ca was varied (since that is the only
%   parameter to affect j_cG_max in an irregular way. Indeed 
%   (1+delta_r)*j_cG_max/B_Ca*F = [j + delta_r*j]/B_Ca*F 
%                               = [j/B_Ca*F] + delta_r*[j/B_Ca*F]
%   So nothing special is needed when varying the normalized j parameter

B_Ca_old = B_Ca/(1+delta_r);
j_cG_max = j_cG_max_normalized*B_Ca_old*F;

j_cG_max_normalized = j_cG_max/(B_Ca*F);






end

