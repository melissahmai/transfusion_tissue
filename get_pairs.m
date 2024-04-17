%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function get_pairs
% Melissa H Mai
%
% Find all cellular pairs in the network (sub-function for build_network.m)
%
% INPUTS
% chi           Connectivity matrix
% tab           Table of cell identities
%
% OUTPUT
% cpairs        List of all pairs in the network
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cpairs = get_pairs(chi,tab)
    cpairs = zeros(sum(chi(:)),2);
    k = 1;
    for i = 1:length(tab.Label)
        if ismember(tab.Label{i}, {'bp','bx','ms','air'})
            continue
        end
        ipair = find(chi(i,:)==1)';
        M = length(ipair);
        cpairs(k:(k+M-1),:) = [repelem(i,M,1) ipair];
        k = k+M;
    end
end