%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script editchi
% Melissa H Mai
%
% Subscript to edit intercellular connection designations for
% build_network.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax = strcmp('ax',tab.Label);
tt = strcmp('tt',tab.Label);
ap = strcmp('ap',tab.Label);
tp = strcmp('tp',tab.Label);
bs = strcmp('bs',tab.Label);
ms = strcmp('ms',tab.Label);
air = strcmp('air',tab.Label);

if ~exist('chi','var') || sum(ismember(tab.Label, {'ax','tt','ap','tp','bs'})) ~= size(chi,1)
    chi = zeros(length(tab.Label));
    % autoconnect bundlesheath
    chi(bs,bs) = [zeros(sum(bs),1) [eye(sum(bs)-1); zeros(1,sum(bs)-1)]];
    
    % start autoconnecting tt
    tt_rows = tab(tt,:);
    for i = 1:length(tt_rows.x)
        r2 = (tt_rows.x-tt_rows.x(i)).^2 + (tt_rows.y-tt_rows.y(i)).^2;
        closest = sort(mink(r2,3),'descend');
        closest = closest(1:2);
        closest_ind = tt_rows.ID(ismember(r2,closest));
        chi(tt_rows.ID(i), closest_ind) = 1;
    end
    
    chi = triu(sign(chi + chi')*1);
end
clf
shownetwork(img,tab,cmap,chi);

loop = true;

fig.Name = 'Select connections';

hold on
act = 'n';
redraw = false;
while strcmpi(act,'n')
    loop = true;
    while loop
        % Select points
        [x,y] = getpts(fig);
        
        if isempty(x) || length(x) == 1
            loop = false;
            break
        end
        
        % Find nearest points
        r2 = (x'-tab.x).^2 + (y'-tab.y).^2;
        pts = zeros(length(x),1);
        for i = 1:length(x)
            pts(i) = find(r2(:,i)==min(r2(:,i)));
            
            if i > 1
                % Register connection in chi
                ind1 = min(pts(1),pts(i));
                ind2 = max(pts(1),pts(i));
                chi(ind1,ind2) = ~logical(chi(ind1,ind2))*1;
                if tt(pts(1)) || ax(pts(1))
                    color = 'blue';
                elseif tt(pts(i)) || ax(pts(i))
                    color = 'blue';
                else
                    color = 'red';
                end
                plot(tab.x([pts(1) pts(i)]), tab.y([pts(1) pts(i)]), '-', 'Color', color)
                if chi(ind1,ind2)==0
                    redraw = true;
                end
            end
        end
        
        if redraw
            shownetwork(img,tab,cmap,chi)
        end
        redraw = false;
    end
    
    fig.Name = groupname;
    shownetwork(img,tab,cmap,chi)
    act = input('Would you like to edit points/network/continue? (p/n/c) ', 's');
    if strcmpi(act,'p')
        editpoints
        act = 'n';
    end
    
end