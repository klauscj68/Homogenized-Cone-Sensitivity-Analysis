function [] = surf_plot(data)
%Plot x,y,z data in 3 space
%   data is a n_pts x 3 array with columns being x,y,z coordinates resp

%% Setup 
x = data(:,1);
y = data(:,2);
z = data(:,3);

%% Make a delaunay triangulation of x-y points to define patches
dt = delaunay(x,y);
n_tri = size(dt,1);

%% Plot surface
F = figure;
hold on

for i=1:n_tri
    fill3(x(dt(i,:)),y(dt(i,:)),z(dt(i,:)),z(dt(i,:)),...
          'EdgeColor','none','FaceColor','interp');
end
colorbar;
view(3);


end

