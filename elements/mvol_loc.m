function [Mvol] = mvol_loc(pts, prism,...
                           R_b,R_t,H)
%Compute \int\phi_i\phi_jdV of shape functions over Eul cone element  
%   Shape functions are tensor products of type affine(x,y) \x affine(z),
%    when pulled back to Lag coordinates cyl elem, which interpolate
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
Mvol = [-(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((0.3e1 / 0.5e1 * z1 ^ 2 + (-0.3e1 / 0.2e1 * H + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - z2 * H / 0.2e1 + z2 ^ 2 / 0.10e2) * R_b ^ 2 + 0.3e1 / 0.2e1 * (-0.4e1 / 0.5e1 * z1 ^ 2 + (H - 0.2e1 / 0.5e1 * z2) * z1 + (H - 0.2e1 / 0.5e1 * z2) * z2 / 0.3e1) * R_t * R_b + 0.3e1 / 0.5e1 * (z1 ^ 2 + z1 * z2 / 0.2e1 + z2 ^ 2 / 0.6e1) * R_t ^ 2) / R_t ^ 2 / H ^ 2 / 0.36e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((0.3e1 / 0.5e1 * z1 ^ 2 + (-0.3e1 / 0.2e1 * H + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - z2 * H / 0.2e1 + z2 ^ 2 / 0.10e2) * R_b ^ 2 + 0.3e1 / 0.2e1 * (-0.4e1 / 0.5e1 * z1 ^ 2 + (H - 0.2e1 / 0.5e1 * z2) * z1 + (H - 0.2e1 / 0.5e1 * z2) * z2 / 0.3e1) * R_t * R_b + 0.3e1 / 0.5e1 * (z1 ^ 2 + z1 * z2 / 0.2e1 + z2 ^ 2 / 0.6e1) * R_t ^ 2) / R_t ^ 2 / H ^ 2 / 0.72e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((0.3e1 / 0.5e1 * z1 ^ 2 + (-0.3e1 / 0.2e1 * H + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - z2 * H / 0.2e1 + z2 ^ 2 / 0.10e2) * R_b ^ 2 + 0.3e1 / 0.2e1 * (-0.4e1 / 0.5e1 * z1 ^ 2 + (H - 0.2e1 / 0.5e1 * z2) * z1 + (H - 0.2e1 / 0.5e1 * z2) * z2 / 0.3e1) * R_t * R_b + 0.3e1 / 0.5e1 * (z1 ^ 2 + z1 * z2 / 0.2e1 + z2 ^ 2 / 0.6e1) * R_t ^ 2) / R_t ^ 2 / H ^ 2 / 0.72e2 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.72e2 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3; -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((0.3e1 / 0.5e1 * z1 ^ 2 + (-0.3e1 / 0.2e1 * H + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - z2 * H / 0.2e1 + z2 ^ 2 / 0.10e2) * R_b ^ 2 + 0.3e1 / 0.2e1 * (-0.4e1 / 0.5e1 * z1 ^ 2 + (H - 0.2e1 / 0.5e1 * z2) * z1 + (H - 0.2e1 / 0.5e1 * z2) * z2 / 0.3e1) * R_t * R_b + 0.3e1 / 0.5e1 * (z1 ^ 2 + z1 * z2 / 0.2e1 + z2 ^ 2 / 0.6e1) * R_t ^ 2) / R_t ^ 2 / H ^ 2 / 0.72e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((0.3e1 / 0.5e1 * z1 ^ 2 + (-0.3e1 / 0.2e1 * H + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - z2 * H / 0.2e1 + z2 ^ 2 / 0.10e2) * R_b ^ 2 + 0.3e1 / 0.2e1 * (-0.4e1 / 0.5e1 * z1 ^ 2 + (H - 0.2e1 / 0.5e1 * z2) * z1 + (H - 0.2e1 / 0.5e1 * z2) * z2 / 0.3e1) * R_t * R_b + 0.3e1 / 0.5e1 * (z1 ^ 2 + z1 * z2 / 0.2e1 + z2 ^ 2 / 0.6e1) * R_t ^ 2) / R_t ^ 2 / H ^ 2 / 0.36e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((0.3e1 / 0.5e1 * z1 ^ 2 + (-0.3e1 / 0.2e1 * H + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - z2 * H / 0.2e1 + z2 ^ 2 / 0.10e2) * R_b ^ 2 + 0.3e1 / 0.2e1 * (-0.4e1 / 0.5e1 * z1 ^ 2 + (H - 0.2e1 / 0.5e1 * z2) * z1 + (H - 0.2e1 / 0.5e1 * z2) * z2 / 0.3e1) * R_t * R_b + 0.3e1 / 0.5e1 * (z1 ^ 2 + z1 * z2 / 0.2e1 + z2 ^ 2 / 0.6e1) * R_t ^ 2) / R_t ^ 2 / H ^ 2 / 0.72e2 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.72e2 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3; -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((0.3e1 / 0.5e1 * z1 ^ 2 + (-0.3e1 / 0.2e1 * H + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - z2 * H / 0.2e1 + z2 ^ 2 / 0.10e2) * R_b ^ 2 + 0.3e1 / 0.2e1 * (-0.4e1 / 0.5e1 * z1 ^ 2 + (H - 0.2e1 / 0.5e1 * z2) * z1 + (H - 0.2e1 / 0.5e1 * z2) * z2 / 0.3e1) * R_t * R_b + 0.3e1 / 0.5e1 * (z1 ^ 2 + z1 * z2 / 0.2e1 + z2 ^ 2 / 0.6e1) * R_t ^ 2) / R_t ^ 2 / H ^ 2 / 0.72e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((0.3e1 / 0.5e1 * z1 ^ 2 + (-0.3e1 / 0.2e1 * H + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - z2 * H / 0.2e1 + z2 ^ 2 / 0.10e2) * R_b ^ 2 + 0.3e1 / 0.2e1 * (-0.4e1 / 0.5e1 * z1 ^ 2 + (H - 0.2e1 / 0.5e1 * z2) * z1 + (H - 0.2e1 / 0.5e1 * z2) * z2 / 0.3e1) * R_t * R_b + 0.3e1 / 0.5e1 * (z1 ^ 2 + z1 * z2 / 0.2e1 + z2 ^ 2 / 0.6e1) * R_t ^ 2) / R_t ^ 2 / H ^ 2 / 0.72e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((0.3e1 / 0.5e1 * z1 ^ 2 + (-0.3e1 / 0.2e1 * H + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - z2 * H / 0.2e1 + z2 ^ 2 / 0.10e2) * R_b ^ 2 + 0.3e1 / 0.2e1 * (-0.4e1 / 0.5e1 * z1 ^ 2 + (H - 0.2e1 / 0.5e1 * z2) * z1 + (H - 0.2e1 / 0.5e1 * z2) * z2 / 0.3e1) * R_t * R_b + 0.3e1 / 0.5e1 * (z1 ^ 2 + z1 * z2 / 0.2e1 + z2 ^ 2 / 0.6e1) * R_t ^ 2) / R_t ^ 2 / H ^ 2 / 0.36e2 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.72e2; -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.72e2 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((z1 ^ 2 / 0.10e2 + (-H / 0.2e1 + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - 0.3e1 / 0.2e1 * z2 * H + 0.3e1 / 0.5e1 * z2 ^ 2) * R_b ^ 2 + R_t * (-0.2e1 / 0.5e1 * z1 ^ 2 + (H - 0.6e1 / 0.5e1 * z2) * z1 + 0.3e1 * (H - 0.4e1 / 0.5e1 * z2) * z2) * R_b / 0.2e1 + R_t ^ 2 * (z1 ^ 2 + 0.3e1 * z1 * z2 + 0.6e1 * z2 ^ 2) / 0.10e2) / R_t ^ 2 / H ^ 2 / 0.36e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((z1 ^ 2 / 0.10e2 + (-H / 0.2e1 + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - 0.3e1 / 0.2e1 * z2 * H + 0.3e1 / 0.5e1 * z2 ^ 2) * R_b ^ 2 + R_t * (-0.2e1 / 0.5e1 * z1 ^ 2 + (H - 0.6e1 / 0.5e1 * z2) * z1 + 0.3e1 * (H - 0.4e1 / 0.5e1 * z2) * z2) * R_b / 0.2e1 + R_t ^ 2 * (z1 ^ 2 + 0.3e1 * z1 * z2 + 0.6e1 * z2 ^ 2) / 0.10e2) / R_t ^ 2 / H ^ 2 / 0.72e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((z1 ^ 2 / 0.10e2 + (-H / 0.2e1 + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - 0.3e1 / 0.2e1 * z2 * H + 0.3e1 / 0.5e1 * z2 ^ 2) * R_b ^ 2 + R_t * (-0.2e1 / 0.5e1 * z1 ^ 2 + (H - 0.6e1 / 0.5e1 * z2) * z1 + 0.3e1 * (H - 0.4e1 / 0.5e1 * z2) * z2) * R_b / 0.2e1 + R_t ^ 2 * (z1 ^ 2 + 0.3e1 * z1 * z2 + 0.6e1 * z2 ^ 2) / 0.10e2) / R_t ^ 2 / H ^ 2 / 0.72e2; -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.72e2 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((z1 ^ 2 / 0.10e2 + (-H / 0.2e1 + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - 0.3e1 / 0.2e1 * z2 * H + 0.3e1 / 0.5e1 * z2 ^ 2) * R_b ^ 2 + R_t * (-0.2e1 / 0.5e1 * z1 ^ 2 + (H - 0.6e1 / 0.5e1 * z2) * z1 + 0.3e1 * (H - 0.4e1 / 0.5e1 * z2) * z2) * R_b / 0.2e1 + R_t ^ 2 * (z1 ^ 2 + 0.3e1 * z1 * z2 + 0.6e1 * z2 ^ 2) / 0.10e2) / R_t ^ 2 / H ^ 2 / 0.72e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((z1 ^ 2 / 0.10e2 + (-H / 0.2e1 + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - 0.3e1 / 0.2e1 * z2 * H + 0.3e1 / 0.5e1 * z2 ^ 2) * R_b ^ 2 + R_t * (-0.2e1 / 0.5e1 * z1 ^ 2 + (H - 0.6e1 / 0.5e1 * z2) * z1 + 0.3e1 * (H - 0.4e1 / 0.5e1 * z2) * z2) * R_b / 0.2e1 + R_t ^ 2 * (z1 ^ 2 + 0.3e1 * z1 * z2 + 0.6e1 * z2 ^ 2) / 0.10e2) / R_t ^ 2 / H ^ 2 / 0.36e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((z1 ^ 2 / 0.10e2 + (-H / 0.2e1 + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - 0.3e1 / 0.2e1 * z2 * H + 0.3e1 / 0.5e1 * z2 ^ 2) * R_b ^ 2 + R_t * (-0.2e1 / 0.5e1 * z1 ^ 2 + (H - 0.6e1 / 0.5e1 * z2) * z1 + 0.3e1 * (H - 0.4e1 / 0.5e1 * z2) * z2) * R_b / 0.2e1 + R_t ^ 2 * (z1 ^ 2 + 0.3e1 * z1 * z2 + 0.6e1 * z2 ^ 2) / 0.10e2) / R_t ^ 2 / H ^ 2 / 0.72e2; -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.144e3 -((0.3e1 / 0.10e2 * z1 ^ 2 + (-H + 0.2e1 / 0.5e1 * z2) * z1 + H ^ 2 - z2 * H + 0.3e1 / 0.10e2 * z2 ^ 2) * R_b ^ 2 + (-0.3e1 / 0.5e1 * z1 ^ 2 + (H - 0.4e1 / 0.5e1 * z2) * z1 + z2 * (H - 0.3e1 / 0.5e1 * z2)) * R_t * R_b + 0.3e1 / 0.10e2 * (z1 ^ 2 + 0.4e1 / 0.3e1 * z1 * z2 + z2 ^ 2) * R_t ^ 2) * (p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) / R_t ^ 2 / H ^ 2 / 0.72e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((z1 ^ 2 / 0.10e2 + (-H / 0.2e1 + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - 0.3e1 / 0.2e1 * z2 * H + 0.3e1 / 0.5e1 * z2 ^ 2) * R_b ^ 2 + R_t * (-0.2e1 / 0.5e1 * z1 ^ 2 + (H - 0.6e1 / 0.5e1 * z2) * z1 + 0.3e1 * (H - 0.4e1 / 0.5e1 * z2) * z2) * R_b / 0.2e1 + R_t ^ 2 * (z1 ^ 2 + 0.3e1 * z1 * z2 + 0.6e1 * z2 ^ 2) / 0.10e2) / R_t ^ 2 / H ^ 2 / 0.72e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((z1 ^ 2 / 0.10e2 + (-H / 0.2e1 + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - 0.3e1 / 0.2e1 * z2 * H + 0.3e1 / 0.5e1 * z2 ^ 2) * R_b ^ 2 + R_t * (-0.2e1 / 0.5e1 * z1 ^ 2 + (H - 0.6e1 / 0.5e1 * z2) * z1 + 0.3e1 * (H - 0.4e1 / 0.5e1 * z2) * z2) * R_b / 0.2e1 + R_t ^ 2 * (z1 ^ 2 + 0.3e1 * z1 * z2 + 0.6e1 * z2 ^ 2) / 0.10e2) / R_t ^ 2 / H ^ 2 / 0.72e2 -(p11 * (p22 - p23) + (-p21 + p23) * p12 - (-p21 + p22) * p13) * (z1 - z2) * ((z1 ^ 2 / 0.10e2 + (-H / 0.2e1 + 0.3e1 / 0.10e2 * z2) * z1 + H ^ 2 - 0.3e1 / 0.2e1 * z2 * H + 0.3e1 / 0.5e1 * z2 ^ 2) * R_b ^ 2 + R_t * (-0.2e1 / 0.5e1 * z1 ^ 2 + (H - 0.6e1 / 0.5e1 * z2) * z1 + 0.3e1 * (H - 0.4e1 / 0.5e1 * z2) * z2) * R_b / 0.2e1 + R_t ^ 2 * (z1 ^ 2 + 0.3e1 * z1 * z2 + 0.6e1 * z2 ^ 2) / 0.10e2) / R_t ^ 2 / H ^ 2 / 0.36e2;];

%% Check for positivity of Maple's chosen Jacobian
% Maple Jcb was J_F*Jbary*Jzed = ( >= 0)*Jbary*(z2-z1)
%  Sign is determined by Jbary
Jbary = p11 * (p22 - p23) + (-p21 + p23) * p12 + p13 * (p21 - p22);

% If Jcb negative, multiply Mvol by -1 to adjust
if Jbary < 0
    Mvol = -1*Mvol;
end

end
