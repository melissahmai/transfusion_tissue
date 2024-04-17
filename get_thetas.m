%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function get_thetas
% Melissa H Mai
%
% Obtain the polar angle for each cell in the network, with theta > 0
% adaxial and theta < 0 abaxial
%
% INPUT
% tab       Table of node identities and coordinates (from [config]_coords.csv)
%
% OUTPUT
% thetas    Vector of angles for each node
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function thetas = get_thetas(tab)
    % Get x & y coordinates
    coords = [tab.x tab.y];
    
    % Isolate bundle sheath coordinates
    bscoords = coords(strcmp(tab.Label,'bs'),:);

    % Find the center of the domain using the bundle sheath
    bscenter = mean(coords(strcmp(tab.Label,'ax')|strcmp(tab.Label,'ap'),:)) - bscoords(end,:);
    
    % Find the center of the inner edge (axis of symmetry)
    center = mean(bscoords([1 end],:)) - bscoords(end,:);
    center = ((bscenter * center')/(center * center') .* center) + bscoords(end,:);

    % Unit vector 1 along the inner edge
    n1 = bscoords(end,:)-bscoords(1,:);
    n1 = n1/norm(n1);
    
    % Unit vector 2 just orthogonal to that
    n2 = [n1(2) -n1(1)];
    
    % Unit vector 3 from edge to flank (approx) through center
    n3 = bscenter + bscoords(end,:) - center;
    n3 = n3/norm(n3);
    
    % If n2 and n3 are in opposite directions, flip n2
    if sign(n3(1))~=sign(n2(1))
        n2 = -n2;
    end

    % Find xylem and phloem coordinates relative to edge center
    bxcoords = coords(strcmp(tab.Label,'ax'),:)-center;
    bxcoords = bxcoords/norm(bxcoords);
    bpcoords = coords(strcmp(tab.Label,'ap'),:)-center;
    bpcoords = bpcoords/norm(bpcoords);

    % Get preliminary orientation from dot product with n1 (inner edge)
    bptheta = acos(bpcoords*n1');
    bxtheta = acos(bxcoords*n1');

    % Flip n1 if necessary so that phloem is abaxial (negative)
    if bxtheta-bptheta > 0
        n1 = -n1;
    end

    dr = (coords - center);
    dr = dr./(sum(dr.^2,2)).^0.5;
    thetas = real(asin(dr*n1'));
end