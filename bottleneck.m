%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function bottleneck
% Melissa H Mai
%
% Build minimal transfusion tissue network models.
%
% INPUTS (mandatory)
% B                 Number of bundle sheath cells
% L                 Number of layers of transfusion tissue
% F                 Funneling factor, which describes how many cells are removed
%                       with each interior layer of transfusion tissue
%
%
% INPUTS (optional)
% 'writetofile',    Boolean value whether coordinates and connectivity
%   'write', 'save'     information are written to files. Filenames will be
%                       b[B]l[L]f[F]_coords.csv and b[B]l[L]f[F]_connections.csv, 
%                       with [B], [L], [F] defined above.
%                       Default: false
%
% 'showfigure',     Boolean value whether to show a diagram of the
%   'showfig'           resulting network
%                       Default: true
%
% 'axap',           Boolean value whether to connect ax and ap
%   'connectvasc'       Default: true
%
% 'markersize'      Marker size for the shown diagram if 'showfig' is true
%                       Default: 15
%
%
% OUTPUTS
% chi               Connectivity matrix
% IN                Identity list
% cx                x-coordinates
% cy                y-coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Building bottlenecking systems
function [chi,IN,cx,cy] = bottleneck(B,L,F,varargin)
    load cmap
    for i = 1:2:length(varargin)
        this = varargin{i+1};
        switch lower(varargin{i})
            case {'writetofile', 'write', 'save'}
                writetofile = this;
            case {'showfig','showfigure'}
                showfig = this;
            case {'markersize'}
                markersize = this;
            case {'axap', 'connectvasc'}
                axap = this;
        end
    end
    if ~exist('writetofile','var'); writetofile = false; end
    if ~exist('showfig','var'); showfig = true; end
    if ~exist('markersize','var'); markersize = 15; end
    if ~exist('axap','var'); axap = true; end
    
    Lm = ceil(B/F);
    edge_exclude = floor(F/2);
    Emod = mod(F,2);
    if L > Lm
        ntp = L*B + 0.5*F*Lm*(1-Lm);
    else
        ntp = L*B - 0.5*F*L*(L-1);
    end

    IN = cellstr(vertcat('ax',repmat('tt',ntp,1),...
        'ap',repmat('tp',ntp,1),...
        repmat('bs',B,1)));

    N = length(IN);
    chi = zeros(N);

    cx = zeros(length(chi),1);
    cx(strcmp(IN,'ax')) = -1.25;
    cy = cx;
    cy(strcmp(IN,'ax')) = -2.5;
    cx(strcmp(IN,'bs')) = 4*(L+1)*ones(B,1);
    cy(strcmp(IN,'bs')) = 4*((0:(B-1))-(B-1)/2)+cy(1)/2;

    ap = find(strcmp(IN,'ap'));

    % Bundle sheath layer
    if B>1
        chi((N-B+1):(end-1), (N-B+2):end) = eye(B-1);
    end
    prev = 1;
    if L == 0
        to_av = unique([floor((B+1)/2) ceil((B+1)/2)]);
        chi(1,to_av+2) = 1;
        chi(2,to_av+2) = 1;
    end
    for l = 1:L
        % #cells in layer l
        if l>Lm
            ncl = B;
        else
            ncl = B - (min(L,Lm) - l)*F;
        end

        if ncl > 0
            % Tracheid
            chi((1:(ncl-1))+prev,(2:ncl)+prev) = eye(ncl-1);

            % Parenchyma
            chi((1:(ncl-1))+prev+ap-1,(2:ncl)+prev+ap-1) = eye(ncl-1);
        end

        % Connect tracheids and parenchyma within layer
        chi((1:ncl)+prev, (1:ncl)+prev+ap-1) = chi((1:ncl)+prev, (1:ncl)+prev+ap-1)+eye(ncl);

        % Connect to previous layer
        if l > Lm
            chi((1:ncl)+prev-ncl,(1:ncl)+prev) = chi((1:ncl)+prev-ncl,(1:ncl)+prev)+eye(ncl);
            chi((1:ncl)+prev-ncl+ap-1,(1:ncl)+prev+ap-1) = chi((1:ncl)+prev-ncl+ap-1,(1:ncl)+prev+ap-1)+eye(ncl);
        elseif l > 1 && ncl >= 2
            if Emod == 1
                cells = ((edge_exclude+1):(ncl-edge_exclude-1))+prev;
                chi(cells-ncl+1+edge_exclude, cells) = chi(cells-ncl+1+edge_exclude, cells)+eye(ncl-F);
                chi(cells-ncl+1+edge_exclude, cells+1) = chi(cells-ncl+1+edge_exclude, cells+1)+eye(ncl-F);
                chi(cells-ncl+ap+edge_exclude, cells+ap-1) = chi(cells-ncl+ap+edge_exclude, cells+ap-1)+eye(ncl-F);
                chi(cells-ncl+ap+edge_exclude, cells+ap) = chi(cells-ncl+ap+edge_exclude, cells+ap)+eye(ncl-F);
            else
                cells = ((edge_exclude+1):(ncl-edge_exclude))+prev;
                chi(cells-ncl+edge_exclude, cells) = chi(cells-ncl+edge_exclude, cells)+eye(ncl-F);
                chi(cells-ncl-1+ap+edge_exclude, cells+ap-1) = chi(cells-ncl-1+ap+edge_exclude, cells+ap-1)+eye(ncl-F);
            end
        elseif ncl == 1
            chi(prev,prev+1)=1;
            chi(prev+ap-1,prev+ap)=1;
%             chi(1,2) = 1;
%             to_cent = unique([floor(F/2) ceil(F/2)]);
%             chi(2,to_cent+3) = 1;
%             chi(ap,ap+1) = 1;
%             chi(ap+1,ap+to_cent+2) = 1;
        elseif l == 1
            to_av = unique([floor((ncl+1)/2) ceil((ncl+1)/2)]);
            chi(1,to_av+1) = 1;
            chi(ap,to_av+ap) = 1;
        end


        % Last layer, connect to bundle sheath
        if l == L
            chi((1:ncl)+prev+ap-1,(N-B+1):end) = eye(ncl);
            chi((1:ncl)+prev,(N-B+1):end) = eye(ncl);
        end
        cx((1:ncl)+prev) = 4*l*ones(ncl,1) + cx(1);
        cy((1:ncl)+prev) = 4*((0:(ncl-1)) - (ncl-1)/2) + cy(1);

        cx((1:ncl)+prev+ap-1) = 4*l*ones(ncl,1);
        cy((1:ncl)+prev+ap-1) = 4*((0:(ncl-1)) - (ncl-1)/2);

        
        prev = prev + ncl;
    end
    
    % If connecting vasculature
    if axap
        chi(1,strcmp(IN,'ap')) = 1;
    end
    
    chi = sign(chi)*1;
    
    if writetofile
        ncol = max(sum(chi,2));
        fid = fopen(sprintf('output/b%dl%df%d_connections.csv',B,L,F), 'w');
        
        for i = 1:length(IN)
            fprintf(fid,'%s',IN{i});
            fprintf(fid,',%d',find(chi(i,:)==1));
            fprintf(fid,repelem(',',ncol-sum(chi(i,:))));
            if i < length(IN)
                fprintf(fid,'\n');
            end
        end
        fclose(fid);
        tab = table(IN,cx,cy,(1:N)','VariableNames',{'Var1','Var2','Var3','Var4'});
        mstab = table(cellstr(repmat('ms',B,1)), 4*(L+2)*ones(B,1), 4*((0:(B-1))-(B-1)/2)'+cy(1)/2, ((N+1):(N+B))');
        airtab = table(cellstr(repmat('air',B,1)), 4*(L+3)*ones(B,1), 4*((0:(B-1))-(B-1)/2)'+cy(1)/2, -200-((N-B+1):N)');
        tab = vertcat(tab, cell2table({'bp',-4,0,0;'bx',cx(1)-4,cy(1),-1}), mstab, airtab);
        tab.Properties.VariableNames = {'Label','x','y','ID'};
        writetable(tab, sprintf('output/b%dl%df%d_coords.csv',B,L,F))
        fprintf('Info written to output/b%dl%df%d_coords.csv\n',B,L,F)
    end
    
    if showfig
%         clf
        for i = 1:(N-1)
            for j = i:N
                if chi(i,j)
                    plot(cx([i j]), cy([i j]), '-', 'Color', [1 1 1]*.7)
                    hold on
                end
            end
        end
        % Axial vasculature
        plot(0,0,'o','MarkerFaceColor',cmap.ap,'Color',cmap.ap,'MarkerSize',markersize);
        hold on
        plot(cx(1),cy(1),'o','MarkerFaceColor',cmap.ax,'Color',cmap.ax,'MarkerSize',markersize);
        
        % Bundle sheath
        plot(cx(strcmp(IN,'bs')), cy(strcmp(IN,'bs')), 's', 'Color', cmap.bs, 'MarkerFaceColor', cmap.bs, 'MarkerSize', markersize);
        
        % Transfusion tissue
        plot(cx(strcmp(IN,'tt')), cy(strcmp(IN,'tt')), 'd', 'Color', cmap.tt, 'MarkerFaceColor', cmap.tt, 'MarkerSize', markersize);
        plot(cx(strcmp(IN,'tp')), cy(strcmp(IN,'tp')), 'd', 'Color', cmap.tp, 'MarkerFaceColor', cmap.tp, 'MarkerSize', markersize);
        hold off
        set(gca,'Visible','off')
    end
end