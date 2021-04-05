function [angle] = uccos(p)
%Unit circle cosine.  Given points p, give their [0,2*pi) angles
%   p is a 2 x n_pts array.  Any ||p|| <=  returns angle NaN
%   angle is a number in [0,2*pi) representing its position on unit circle
x = p(1,:);
y = p(2,:);

%% Normalize points to unit circle
radius = sqrt(x.^2 + y.^2);
radius(radius == 0) = NaN;

ux = x./radius;

%% Compute acos of all points
angle = acos(ux);

%% Determine quadrant points that must be shifted
%quad1 = find((x >= 0).*(y >= 0));
%quad2 = find((x < 0).*(y >= 0));
quad3 = find((x <= 0).*(y < 0));
quad4 = find((x > 0).*(y < 0));

%% Adjust angle by appropriate quadrants
% Quadrant 3
angle(quad3) = 2*pi - angle(quad3);

% Quadrant 4
angle(quad4) = 2*pi - angle(quad4);


end

