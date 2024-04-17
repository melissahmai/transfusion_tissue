%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script editpoints
% Melissa H Mai
%
% Subscript to edit cell identifications for
% build_network.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('act','var')
    act = '';
end

if ~exist('chi','var')
    chi = 0;
end
shownetwork(img, tab, cmap, chi)
ax = strcmp('ax',tab.Label);
tt = strcmp('tt',tab.Label);
ap = strcmp('ap',tab.Label);
tp = strcmp('tp',tab.Label);
bs = strcmp('bs',tab.Label);
ms = strcmp('ms',tab.Label);
air = strcmp('air',tab.Label);

while ~strcmpi(act,'c')
    fig.Name = groupname;
    act = input('Would you like to add/delete/move/continue? (a/d/m/c) ', 's');
    
    if strcmpi(act,'a')
        cont = true;
        while cont
            % Ask for cell type to add
            celltype = input('Which cell type? (ax/ap/tt/tp/bs/c) ', 's');
            celltype = lower(celltype);
            
            % Check input
            if ~ismember(celltype, {'ax','ap','tt','tp','bs','c'})
                disp('Invalid input')
            elseif strcmp(celltype,'c')
                cont = false;
            else
                fig.Name = strjoin({'Add ', celltype, ' points'});
                cont2 = true;
                
                while cont2
                    [xadd,yadd] = getpts(fig);
                    
                    % If empty, end
                    if isempty(xadd)
                        cont2 = false;
                        continue
                    end
                    
                    % Find nearest points
                    r2 = (xadd'-tab.x(strcmp(celltype,tab.Label))).^2 + ...
                        (yadd'-tab.y(strcmp(celltype,tab.Label))).^2;
                    
                    for i = 1:length(xadd)
                        nearest = find(ismember(r2(:,i),mink(r2(:,i),2)));
                        newpt = nearest(2) + min(tab.ID(strcmp(celltype,tab.Label))) - 1;
                        
                        N = length(tab.ID);
                        tab((newpt+1):(N+1),:) = tab(newpt:end,:);
                        tab(newpt,:) = table({celltype},xadd(i),yadd(i),i);
                        
                        % Reassign IDs
                        tab.ID = (1:length(tab.ID))';
                        nbs = sum(strcmp(tab.Label,'bs'));
                        N = sum(ismember(tab.Label, {'ax','ap','tt','tp','bs'}));
                        M = sum(strcmp(tab.Label,'ms'));
                        if M > 0
                            tab.ID((N+M+1):end) = -200-N+nbs - (1:nbs)';
                        end
                        
                        % Adjust chi if it exists
                        if exist('chi','var') && length(chi)>1
                            chi((newpt+1):(end+1),:) = chi(newpt:end,:);
                            chi(:,(newpt+1):(end+1)) = chi(:,newpt:end);
                            chi(newpt,:) = 0;
                            chi(:,newpt) = 0;
                        end
                        
                        hold off
                        shownetwork(img, tab, cmap, chi);
                    end
                end
                fig.Name = groupname;
            end
        end
    elseif strcmpi(act,'d')
        cont = true;
        while cont
            fig.Name = 'Delete points';
            % Select points
            [xdel,ydel] = getpts(fig);
            
            % If empty, end
            if isempty(xdel)
                cont = false;
                continue
            end
            
            % Find nearest points
            r2 = (xdel'-tab.x).^2 + (ydel'-tab.y).^2;
            pts = zeros(length(xdel),1);
            for i = 1:length(xdel)
                pts(i) = find(r2(:,i)==min(r2(:,i)),1);
            end
            
            % Remove points
            N = sum(ismember(tab.Label, {'ax','ap','tt','tp','bs'}));
            nbs = sum(strcmp(tab.Label,'bs'));
            tab(ismember(tab.ID,pts),:) = [];
            tab(ismember(tab.ID,pts(pts>N-nbs)+nbs),:) = [];
            tab(ismember(tab.ID,-200-pts(pts>N-nbs)),:) = [];
            
            % Reassign IDs
            tab.ID = (1:length(tab.ID))';
            N = sum(ismember(tab.Label, {'ax','ap','tt','tp','bs'}));
            nbs = sum(strcmp(tab.Label,'bs'));
            tab.ID(strcmp(tab.Label,'bp')) = 0;
            tab.ID(strcmp(tab.Label,'bx')) = -1;
            if sum(strcmp(tab.Label,'air'))>0
                tab.ID(strcmp(tab.Label,'air')) = -200-N+nbs - (1:nbs);
            end
%             
%             M = sum(ismember(tab.Label, {'ms','bp','bx'}));
%             if M > 0
%                 tab.ID((N+M+1):end) = -200-N+nbs - (1:nbs)';
%             end
            
            % Adjust chi if it exists
            if exist('chi','var') && length(chi)>1
                chi(pts,:) = [];
                chi(:,pts) = [];
            end
            
            hold off
            shownetwork(img,tab,cmap, chi);
        end
    elseif strcmpi(act,'m')
        cont = true;
        fig.Name = 'Move points by selecting a point, then its new location, then Enter';
        while cont
            [xmov,ymov] = getpts(fig);
            if isempty(xmov)
                cont = false;
                continue
            end
            
            % Find nearest point
            r2 = (xmov(1)'-tab.x).^2 + (ymov(1)'-tab.y).^2;
            pt = find(r2==min(r2));
            
            tab.x(pt) = xmov(2);
            tab.y(pt) = ymov(2);
            
            hold off
            shownetwork(img, tab, cmap, chi)
            
        end
    end
end