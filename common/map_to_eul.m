function [pts_eul] = map_to_eul(pts,...
                                R_b,R_t,H)
%Maps pts from Lag cylinder domain to Eul cone domain
%   pts is a 3 by n_pts array
%   Remaining args are geometric parameters for cone
%   pts_eul is a 3 by n_pts array

%% Initialize parameters
n_pts = size(pts,2);

%% Compute the scaling coefficient as function of z
lambda = R_b/R_t + (1-R_b/R_t)/H*pts(3,:);

%% Map the pts
M = [lambda;...
    lambda;...
    ones(1,n_pts)];

pts_eul = M.*pts;

end

