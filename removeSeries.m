%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function removeSeries
% Melissa H Mai
%
% Removes a cell type from the get_pathfigs figure.
%
% INPUTS
% thisax        Axis handle from get_pathfigs
% cells         Cell array of cell types to remove
%                   within {'ax', 'ap', 'tt', 'tp', 'bs', 'ms}
%
% OUTPUT
% thisax        Axis handle for the modified figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function thisax = removeSeries(thisax,cells)

    to_delete = [];
    for i = 1:length(thisax.Children)
        this = thisax.Children(i);
        if any(strcmp(this.DisplayName,cells))
            to_delete = [to_delete i];
        end
    end
    delete(thisax.Children(to_delete))
end