function [M] = mvol(pts,prisms,...
                    R_b,R_t,H)
%Build (\int phi_i*phi_j dx) across cone |_|prisms

% Set path to include elements directory
oldpath = addpath('..\elements\');

%% Needed parameters
n_pts = size(pts,2);
n_prisms = size(prisms,2);

%% Initialize lists used to build sparse M
% Number of nodal pairs we will cycle over while redundant (i,j) values
% will be accumulated by sparse
n_pairs = 36*n_prisms;

I = zeros(n_pairs,1);
J = zeros(n_pairs,1);
V = zeros(n_pairs,1);

%% Loop over elements and build lists
for k=1:n_prisms
    % Slice in needed values
    prism = prisms(:,k);
    
    % Compute Mvol_loc
    Mvol_loc = mvol_loc(pts,prism,R_b,R_t,H);
    
    % Write these values in list format
    % Order of indices is compatible across I,J, V because list order is 
    % first down array column then onto next column
    
    %  Array value is i-index there
    ram_I = repmat([1;2;3;4;5;6],1,6);
    ram_I = ram_I(:);
    %  Now switch from local node indexing to mesh indexing
    ram_I = prism(ram_I);
    
    %  Array value is j-index there
    ram_J = repmat([1 2 3 4 5 6],6,1);
    ram_J = ram_J(:);
    %  Now switch from local node indexing to mesh indexing
    ram_J = prism(ram_J);
    
    ram_V = Mvol_loc(:);
    
    % Write these values into I,J,V
    I((k-1)*36+1:k*36) = ram_I;
    J((k-1)*36+1:k*36) = ram_J;
    V((k-1)*36+1:k*36) = ram_V;
end

%% Build M
M = sparse(I,J,V,n_pts,n_pts);

%% Return path to normal
path(oldpath)

end


