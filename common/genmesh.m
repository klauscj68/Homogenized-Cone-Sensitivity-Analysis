function [pts,prism,faces_sl,faces_fo,...
          n_pts,n_prism,n_fsl,n_ffo] = ...
          genmesh(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,sigma,...
                  n_sez,taglia,tol_R,tol_angle)
%Generate mesh lists for Lag 3D cylinder 
%   Detailed explanation goes here

%% Generate 2D cross section mesh

% Define the circle geometry
gd = [1;0;0;R_t];
ns = 'disc';
ns = ns';
sf = 'disc';
disc_obj = decsg(gd,sf,ns);

% Generate initial cross section mesh
%  May view cross section mesh with pdeplot(cross)
[p,e,t] = initmesh(disc_obj,'Hmax',taglia);

cross_pts = p(1:2,:);
cross_tri = t(1:3,:);

n_cpts = size(cross_pts,2);
n_tri = size(cross_tri,2);

%% Call 2D edge and boundary lists
[edge,bdy] = edgelists(cross_pts,cross_tri);
n_edge = size(edge,2);
n_bdy  = size(bdy,2);


%% Refine the mesh near rim and build edge lists
% Initialize refining process
n_refine = 0;
tol_near = .2*R_t;

for i=1:n_refine
% Find the triangles whose midpoints are close to the transition
%  Compute the midpoints
cross_midpts = [1/3*(cross_pts(1,cross_tri(1,:)) + cross_pts(1,cross_tri(2,:)) + cross_pts(1,cross_tri(3,:)));...
                1/3*(cross_pts(2,cross_tri(1,:)) + cross_pts(2,cross_tri(2,:)) + cross_pts(2,cross_tri(3,:)))];
            
%  Test for midpoints close to rim
%   Distance of midpoints to rim
prox = zeros(1,n_tri);
for j=1:n_tri
    prox(:,j) = norm(cross_midpts(:,j));
end

prox_test = prox >= R_t - tol_near;
rtri = find(prox_test);

[p,e,t] = refinemesh(disc_obj,p,e,t,rtri');

% Update the needed variables
cross_pts = p(1:2,:);
cross_tri = t(1:3,:);

n_cpts = size(cross_pts,2);
n_tri = size(cross_tri,2);
                                 
[edge,bdy] = edgelists(cross_pts,cross_tri);
n_edge = size(edge,2);
n_bdy  = size(bdy,2);
            
end

%% List boundary edges on sliver, and complement on folds
flag_sl = false(1,n_bdy);

% x-y coordinates of points on boundary
%  1st two rows are 1st point
%  2nd two rows are 2nd point
ram = [cross_pts(1:2,edge(1,bdy));...
       cross_pts(1:2,edge(2,bdy))];

% Angles of points on boundary
%   Columns match the corresponding edge of bdy
angle = zeros(2,n_bdy);
angle(1,:) = uccos(ram(1:2,:));
angle(2,:) = uccos(ram(3:4,:));

% Find bdy edges with midpoint in sliver
%   Find bdy edges with their midpoint on the sliver (indicated by >= 1 value) 
angle_test = logical((mean(angle) >= theta_in).*...
                     (mean(angle) <= theta_fin));
flag_sl(...
        angle_test...
        ) = true;

% Find edges with midpoint not in sliver
flag_fo = flag_sl;
flag_fo(isnan(angle_test)) = true;
flag_fo = logical(1 - flag_fo);

% Convert logical index to numerical index
bdy_sl = find(flag_sl);
bdy_fo = find(flag_fo);

%% Build 3D node list
assert(n_sez > 1,'n_sez must be greater than 1');

% Duplicate the x-y coordinates for all sections
pts = repmat(cross_pts,1,n_sez);

% Write the z coordinates
dz = H/(n_sez-1);
zvals = index_shift(n_cpts,n_sez);
zvals = zvals - 1; 
zvals = dz*zvals;

% Augment vertex list for 3D
pts = [pts;...
       zvals];

%% Build 3D prism list
% Top three rows are index of bottom triangular face
% Second three rows are index of top triangular face
% Vertices are toured like output of generateMesh
prism = repmat(cross_tri,2,n_sez-1);

shift = index_shift(n_tri,n_sez-1);
shift = repmat(shift,3,1);
prism = [prism(1:3,:) + (shift - 1)*n_cpts;...
         prism(4:6,:) + shift*n_cpts];


%% Build 3D boundary faces lists
% Build the faces on the sliver
edge_sl = edge(1:2,...
                  bdy(bdy_sl));
n_edge_sl = size(edge_sl,2);
faces_sl = repmat(edge_sl,2,n_sez-1);


shift = index_shift(n_edge_sl,n_sez-1);
shift = repmat(shift,2,1);
faces_sl = [faces_sl(1:2,:) + (shift-1)*n_cpts;...
            faces_sl(3:4,:) + shift*n_cpts];
        
% Build the faces on the folds
edge_fo = edge(1:2,...
                  bdy(bdy_fo));
n_edge_fo = size(edge_fo,2);
faces_fo = repmat(edge_fo,2,n_sez-1);

shift = index_shift(n_edge_fo,n_sez-1);
shift = repmat(shift,2,1);
faces_fo = [faces_fo(1:2,:) + (shift-1)*n_cpts;...
            faces_fo(3:4,:) + shift*n_cpts];

%% Output numbers of 3D pts, prisms, and bdy faces
n_pts = size(pts,2);
n_prism = size(prism,2);
n_fsl = size(faces_sl,2);
n_ffo = size(faces_fo,2);

end

