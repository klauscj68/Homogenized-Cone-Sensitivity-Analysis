function [edge,bdy] = edgelists(pts,tri)
%Given a 2D mesh output the edges list and indices of edges on boundary
%   Detailed explanation goes here
n_tri = size(tri,2);

%% Generate the edge list
% Initialize edge list
%  Top two columns are vertex indices defining edge ordered least -> great
%  Next two rows indicate which triangles an edge belongs.  An edge is on
%  boundary iff it belongs to exactly one triangle
edge = zeros(4,3*n_tri);
n_edge = 0;

% Loop over all triangles and add in new edges
for i=1:n_tri
    % Extract i^th triangle
    %  Note: bc we sort T and write cands respecting that order
    %        edges vertices are automatically sorted too as we write
    %        in new edges usings cands
    T = sort(tri(:,i));
    
    % Locate any previously found edges with both vertices lying in T
    %  dupl will be column indices of edges with both ends in T
    if n_edge~=0
        dupl = ismember(edge(1:2,1:n_edge),T);
        dupl = logical(dupl(1,1:n_edge).*dupl(2,1:n_edge));
        dupl = find(dupl);
    else 
        dupl = [];
    end
    
    
    % Possible new edges contributed by triangle
    cands = [ [T(1);T(2)] [T(1);T(3)] [T(2);T(3)] ];
    
    if ~isempty(dupl)
        n_dupl = size(dupl,2);
        
        % flag_new set to false if ever cand edge is listed
        flag_new = true(1,3);
        for j=1:n_dupl
            edgeInT = dupl(j);
            v1 = edge(1,edgeInT);
            v2 = edge(2,edgeInT);
            
            if (v1 == cands(1,1))&&(v2 == cands(2,1))
                flag_new(1) = false;
                % Record this edge also belongs to triangle i
                edge(4,edgeInT) = i;
            end
            
            if (v1 == cands(1,2))&&(v2 == cands(2,2))
                flag_new(2) = false;
                % Record this edge also belongs to triangle i
                edge(4,edgeInT) = i;
            end
            
            if (v1 == cands(1,3))&&(v2 == cands(2,3))
                flag_new(3) = false;
                % Record this edge also belongs to triangle i
                edge(4,edgeInT) = i;
            end
            
        end
        
        % Update edge lists to include any new edges
        new = find(flag_new);
        n_new = size(new,2);
        for j=1:n_new
            edge(1:3,n_edge+j) = [cands(1,new(j));...
                                  cands(2,new(j));...
                                  i];
        end
        
        % Update count of total edges
        n_edge = n_edge + n_new;
       
    else
        % There were no duplicates and will write in all edges
        edge(1:3,n_edge+1:n_edge+3) = [cands;...
                                       i*ones(1,3)];
                                   
        % Update count of total edges
        n_edge = n_edge + 3;
        
    end
end

% Cut out padded zeros
edge = edge(:,1:n_edge);

%% Generate the bdy list
bdy = find(edge(4,:)==0);

end

