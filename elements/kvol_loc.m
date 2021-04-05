function [Kvol] = kvol_loc(pts, prism,...
                           R_b,R_t,H)
%Compute \intD\phi_iD\phi_jdV of shape functions over Eul cone element  
%   Shape functions are tensor products of type affine(x,y) \x affine(z),
%    when pulled back to Lag coordinates cyl element, which interpolate
%    1 at their associated vertex and 0 at all others.
%   pts is a 3 x n_pts node list over the Lag coordinate cylinder mesh
%   prism is a 6 x 1 array whose first three entries are pts indices for
%    the base and next three are pts indices for the nodes lying directly
%    above the base
%    RMK: intended that prism's nodes are listed in prism as 
%         [x1 x2 x3 x1 x2 x3;...
%          z1 z1 z1 z2 z2 z2];
%         with z1 < z2
%   R_b, R_t, H are the geometric parameters of the cone which our Lag
%   coordinates cylinder maps onto by horizontal scaling

%% Extract the local points for this prism
P = pts(:,prism(:)');

%Reformat for use in Maple exported routine
p11 = P(1,1); p12 = P(1,2); p13 = P(1,3);
p21 = P(2,1); p22 = P(2,2); p23 = P(2,3);

z1 = P(3,1);
z2 = P(3,4);

%% Run Maple generated local assembly
Kvol = [-(z1 - z2) * (p12 ^ 2 - 2 * p12 * p13 + p13 ^ 2 + (p22 - p23) ^ 2) / ((-6 * p21 + 6 * p23) * p12 + (6 * p21 - 6 * p22) * p13 + 6 * p11 * (p22 - p23)) (z1 - z2) * (p13 ^ 2 + (-p11 - p12) * p13 + p23 ^ 2 + (-p21 - p22) * p23 + p12 * p11 + p22 * p21) / ((6 * p21 - 6 * p22) * p13 + (-6 * p11 + 6 * p12) * p23 + 6 * p11 * p22 - 6 * p12 * p21) -(-p12 ^ 2 + (p11 + p13) * p12 - p22 ^ 2 + (p21 + p23) * p22 - p13 * p11 - p23 * p21) * (z1 - z2) / ((-6 * p21 + 6 * p23) * p12 + (6 * p11 - 6 * p13) * p22 - 6 * p11 * p23 + 6 * p13 * p21) -(z1 - z2) * (p12 ^ 2 - 2 * p12 * p13 + p13 ^ 2 + (p22 - p23) ^ 2) / ((-12 * p21 + 12 * p23) * p12 + (12 * p21 - 12 * p22) * p13 + 12 * p11 * (p22 - p23)) (z1 - z2) * (p13 ^ 2 + (-p11 - p12) * p13 + p23 ^ 2 + (-p21 - p22) * p23 + p12 * p11 + p22 * p21) / ((12 * p21 - 12 * p22) * p13 + (-12 * p11 + 12 * p12) * p23 + 12 * p11 * p22 - 12 * p12 * p21) -(-p12 ^ 2 + (p11 + p13) * p12 - p22 ^ 2 + (p21 + p23) * p22 - p13 * p11 - p23 * p21) * (z1 - z2) / ((-12 * p21 + 12 * p23) * p12 + (12 * p11 - 12 * p13) * p22 - 12 * p11 * p23 + 12 * p13 * p21); (z1 - z2) * (p13 ^ 2 + (-p11 - p12) * p13 + p23 ^ 2 + (-p21 - p22) * p23 + p12 * p11 + p22 * p21) / ((6 * p21 - 6 * p22) * p13 + (-6 * p11 + 6 * p12) * p23 + 6 * p11 * p22 - 6 * p12 * p21) -(p11 ^ 2 - 2 * p13 * p11 + p13 ^ 2 + (p21 - p23) ^ 2) * (z1 - z2) / ((6 * p22 - 6 * p23) * p11 + (6 * p21 - 6 * p22) * p13 - 6 * p12 * (p21 - p23)) (z1 - z2) * (p11 ^ 2 + (-p12 - p13) * p11 + p21 ^ 2 + (-p22 - p23) * p21 + p12 * p13 + p22 * p23) / ((6 * p22 - 6 * p23) * p11 + (-6 * p12 + 6 * p13) * p21 + 6 * p12 * p23 - 6 * p13 * p22) (z1 - z2) * (p13 ^ 2 + (-p11 - p12) * p13 + p23 ^ 2 + (-p21 - p22) * p23 + p12 * p11 + p22 * p21) / ((12 * p21 - 12 * p22) * p13 + (-12 * p11 + 12 * p12) * p23 + 12 * p11 * p22 - 12 * p12 * p21) -(p11 ^ 2 - 2 * p13 * p11 + p13 ^ 2 + (p21 - p23) ^ 2) * (z1 - z2) / ((12 * p22 - 12 * p23) * p11 + (12 * p21 - 12 * p22) * p13 - 12 * p12 * (p21 - p23)) (z1 - z2) * (p11 ^ 2 + (-p12 - p13) * p11 + p21 ^ 2 + (-p22 - p23) * p21 + p12 * p13 + p22 * p23) / ((12 * p22 - 12 * p23) * p11 + (-12 * p12 + 12 * p13) * p21 + 12 * p12 * p23 - 12 * p13 * p22); -(-p12 ^ 2 + (p11 + p13) * p12 - p22 ^ 2 + (p21 + p23) * p22 - p13 * p11 - p23 * p21) * (z1 - z2) / ((-6 * p21 + 6 * p23) * p12 + (6 * p11 - 6 * p13) * p22 - 6 * p11 * p23 + 6 * p13 * p21) (z1 - z2) * (p11 ^ 2 + (-p12 - p13) * p11 + p21 ^ 2 + (-p22 - p23) * p21 + p12 * p13 + p22 * p23) / ((6 * p22 - 6 * p23) * p11 + (-6 * p12 + 6 * p13) * p21 + 6 * p12 * p23 - 6 * p13 * p22) -(z1 - z2) * (p11 ^ 2 - 2 * p12 * p11 + p12 ^ 2 + (p21 - p22) ^ 2) / ((6 * p22 - 6 * p23) * p11 + (-6 * p21 + 6 * p23) * p12 + 6 * p13 * (p21 - p22)) -(-p12 ^ 2 + (p11 + p13) * p12 - p22 ^ 2 + (p21 + p23) * p22 - p13 * p11 - p23 * p21) * (z1 - z2) / ((-12 * p21 + 12 * p23) * p12 + (12 * p11 - 12 * p13) * p22 - 12 * p11 * p23 + 12 * p13 * p21) (z1 - z2) * (p11 ^ 2 + (-p12 - p13) * p11 + p21 ^ 2 + (-p22 - p23) * p21 + p12 * p13 + p22 * p23) / ((12 * p22 - 12 * p23) * p11 + (-12 * p12 + 12 * p13) * p21 + 12 * p12 * p23 - 12 * p13 * p22) -(z1 - z2) * (p11 ^ 2 - 2 * p12 * p11 + p12 ^ 2 + (p21 - p22) ^ 2) / ((12 * p22 - 12 * p23) * p11 + (-12 * p21 + 12 * p23) * p12 + 12 * p13 * (p21 - p22)); -(z1 - z2) * (p12 ^ 2 - 2 * p12 * p13 + p13 ^ 2 + (p22 - p23) ^ 2) / ((-12 * p21 + 12 * p23) * p12 + (12 * p21 - 12 * p22) * p13 + 12 * p11 * (p22 - p23)) (z1 - z2) * (p13 ^ 2 + (-p11 - p12) * p13 + p23 ^ 2 + (-p21 - p22) * p23 + p12 * p11 + p22 * p21) / ((12 * p21 - 12 * p22) * p13 + (-12 * p11 + 12 * p12) * p23 + 12 * p11 * p22 - 12 * p12 * p21) -(-p12 ^ 2 + (p11 + p13) * p12 - p22 ^ 2 + (p21 + p23) * p22 - p13 * p11 - p23 * p21) * (z1 - z2) / ((-12 * p21 + 12 * p23) * p12 + (12 * p11 - 12 * p13) * p22 - 12 * p11 * p23 + 12 * p13 * p21) -(z1 - z2) * (p12 ^ 2 - 2 * p12 * p13 + p13 ^ 2 + (p22 - p23) ^ 2) / ((-6 * p21 + 6 * p23) * p12 + (6 * p21 - 6 * p22) * p13 + 6 * p11 * (p22 - p23)) (z1 - z2) * (p13 ^ 2 + (-p11 - p12) * p13 + p23 ^ 2 + (-p21 - p22) * p23 + p12 * p11 + p22 * p21) / ((6 * p21 - 6 * p22) * p13 + (-6 * p11 + 6 * p12) * p23 + 6 * p11 * p22 - 6 * p12 * p21) -(-p12 ^ 2 + (p11 + p13) * p12 - p22 ^ 2 + (p21 + p23) * p22 - p13 * p11 - p23 * p21) * (z1 - z2) / ((-6 * p21 + 6 * p23) * p12 + (6 * p11 - 6 * p13) * p22 - 6 * p11 * p23 + 6 * p13 * p21); (z1 - z2) * (p13 ^ 2 + (-p11 - p12) * p13 + p23 ^ 2 + (-p21 - p22) * p23 + p12 * p11 + p22 * p21) / ((12 * p21 - 12 * p22) * p13 + (-12 * p11 + 12 * p12) * p23 + 12 * p11 * p22 - 12 * p12 * p21) -(p11 ^ 2 - 2 * p13 * p11 + p13 ^ 2 + (p21 - p23) ^ 2) * (z1 - z2) / ((12 * p22 - 12 * p23) * p11 + (12 * p21 - 12 * p22) * p13 - 12 * p12 * (p21 - p23)) (z1 - z2) * (p11 ^ 2 + (-p12 - p13) * p11 + p21 ^ 2 + (-p22 - p23) * p21 + p12 * p13 + p22 * p23) / ((12 * p22 - 12 * p23) * p11 + (-12 * p12 + 12 * p13) * p21 + 12 * p12 * p23 - 12 * p13 * p22) (z1 - z2) * (p13 ^ 2 + (-p11 - p12) * p13 + p23 ^ 2 + (-p21 - p22) * p23 + p12 * p11 + p22 * p21) / ((6 * p21 - 6 * p22) * p13 + (-6 * p11 + 6 * p12) * p23 + 6 * p11 * p22 - 6 * p12 * p21) -(p11 ^ 2 - 2 * p13 * p11 + p13 ^ 2 + (p21 - p23) ^ 2) * (z1 - z2) / ((6 * p22 - 6 * p23) * p11 + (6 * p21 - 6 * p22) * p13 - 6 * p12 * (p21 - p23)) (z1 - z2) * (p11 ^ 2 + (-p12 - p13) * p11 + p21 ^ 2 + (-p22 - p23) * p21 + p12 * p13 + p22 * p23) / ((6 * p22 - 6 * p23) * p11 + (-6 * p12 + 6 * p13) * p21 + 6 * p12 * p23 - 6 * p13 * p22); -(-p12 ^ 2 + (p11 + p13) * p12 - p22 ^ 2 + (p21 + p23) * p22 - p13 * p11 - p23 * p21) * (z1 - z2) / ((-12 * p21 + 12 * p23) * p12 + (12 * p11 - 12 * p13) * p22 - 12 * p11 * p23 + 12 * p13 * p21) (z1 - z2) * (p11 ^ 2 + (-p12 - p13) * p11 + p21 ^ 2 + (-p22 - p23) * p21 + p12 * p13 + p22 * p23) / ((12 * p22 - 12 * p23) * p11 + (-12 * p12 + 12 * p13) * p21 + 12 * p12 * p23 - 12 * p13 * p22) -(z1 - z2) * (p11 ^ 2 - 2 * p12 * p11 + p12 ^ 2 + (p21 - p22) ^ 2) / ((12 * p22 - 12 * p23) * p11 + (-12 * p21 + 12 * p23) * p12 + 12 * p13 * (p21 - p22)) -(-p12 ^ 2 + (p11 + p13) * p12 - p22 ^ 2 + (p21 + p23) * p22 - p13 * p11 - p23 * p21) * (z1 - z2) / ((-6 * p21 + 6 * p23) * p12 + (6 * p11 - 6 * p13) * p22 - 6 * p11 * p23 + 6 * p13 * p21) (z1 - z2) * (p11 ^ 2 + (-p12 - p13) * p11 + p21 ^ 2 + (-p22 - p23) * p21 + p12 * p13 + p22 * p23) / ((6 * p22 - 6 * p23) * p11 + (-6 * p12 + 6 * p13) * p21 + 6 * p12 * p23 - 6 * p13 * p22) -(z1 - z2) * (p11 ^ 2 - 2 * p12 * p11 + p12 ^ 2 + (p21 - p22) ^ 2) / ((6 * p22 - 6 * p23) * p11 + (-6 * p21 + 6 * p23) * p12 + 6 * p13 * (p21 - p22));];

%% Check for positivity of Maple's chosen Jacobian
% Maple Jcb was J_F*Jbary*Jzed = ( >= 0)*Jbary*(z2-z1)
%  Sign is determined by Jbary
Jbary = p11 * (p22 - p23) + (-p21 + p23) * p12 + p13 * (p21 - p22);

% If Jcb negative, multiply Mvol by -1 to adjust
if Jbary < 0
    Kvol = -1*Kvol;
end

end

