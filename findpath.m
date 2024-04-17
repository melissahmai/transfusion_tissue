%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function findpath
% Melissa H Mai
%
% Find shortest path from a start point to each cell in the network
%
% INPUTS
% chi           Connectivity matrix
% start         Index for the reference point
%
% OUTPUT
% distances     Vector showing the length of the shortest path from each point
%                   to the reference (reference to reference = 0).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function distances = findpath(chi,start)
    % Get fully symmetric chi matrix
    chi = triu(chi) + triu(chi)';
    
    % Initialize chain
    chain = cell(size(chi,1),1);

    % Find the first nodes adjacent to the start node
    chain{1} = find(sum(chi(start(:),:),1)>0);
    for i = 2:length(chain)
        % Identify the next set of adjacent nodes, removing ones that have
        % already been identified
        chain{i} = setdiff(find(sum(chi(chain{i-1},:),1)>0), horzcat(chain{1:(i-1)},start(:)'));
        if (length(horzcat(chain{1:i}))+length(start)) == size(chi,1)
            % Finish once all nodes have been found
            chain = chain(1:i);
            break
        end
    end

    % Convert chain to distances
    distances = zeros(size(chi,1),1);
    for i = 1:length(chain)
        distances(chain{i}) = i;
    end
    distances(start) = 0;
end