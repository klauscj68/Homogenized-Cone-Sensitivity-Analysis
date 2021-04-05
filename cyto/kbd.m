function [Kbd] = kbd(pts,rects,...
                    R_b,R_t,H)
%Build (\int D_LBphi_i*D_LBphi_j dsigma) across dcone |_|rects

% Set path to include elements directory
oldpath = addpath('..\elements\');

%% Needed parameters
n_pts = size(pts,2);
n_rects = size(rects,2);

%% Initialize lists used to build sparse K
% Number of nodal pairs we will cycle over while redundant (i,j) values
% will be accumulated by sparse
n_pairs = 16*n_rects;

I = zeros(n_pairs,1);
J = zeros(n_pairs,1);
V = zeros(n_pairs,1);

%% Loop over elements and build lists
for k=1:n_rects
    % Slice in needed values
    rect = rects(:,k);
    
    % Compute Kbd_loc
    Kbd_loc = kbd_loc(pts,rect,R_b,R_t,H);
    
    % Write these values in list format
    % Order of indices is compatible across I,J, V because list order is 
    % first down array column then onto next column
    
    %  Array value is i-index there
    ram_I = repmat([1;2;3;4],1,4);
    ram_I = ram_I(:);
    %  Now switch from local node indexing to mesh indexing
    ram_I = rect(ram_I);
    
    %  Array value is j-index there
    ram_J = repmat([1 2 3 4],4,1);
    ram_J = ram_J(:);
    %  Now switch from local node indexing to mesh indexing
    ram_J = rect(ram_J);
    
    ram_V = Kbd_loc(:);
    
    % Write these values into I,J,V
    I((k-1)*16+1:k*16) = ram_I;
    J((k-1)*16+1:k*16) = ram_J;
    V((k-1)*16+1:k*16) = ram_V;
end

%% Build Kbd
Kbd = sparse(I,J,V,n_pts,n_pts);

%% Return path to normal
path(oldpath)

end

