%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function findpath_edge
% Melissa H Mai
%
% Find shortest path from a start point to each connection in the network
%
% INPUTS
% chi           Connectivity matrix
% start         Index for the reference point
%
% OUTPUT
% distances     Matrix showing the length of the shortest path from each connection
%                   to the reference (reference to reference = 0).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Find shortest path (edges)
function distances = findpath_edge(chi,start)
    dnode = findpath(chi,start);

    distances = chi*0 - 1;
    
    distances(dnode==0,:) = 0;

    for i = 1:max(dnode)
        distances(dnode==i,:) = i*chi(dnode==i,:);
    end
    
    upper = triu(distances);
    lower = tril(distances)';
    distances = min(upper,lower);
end