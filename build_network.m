%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script build_network
% Melissa H Mai
%
% Builds transfusion tissue networks for tt_model. Loads in cross-sectional
% images of the needle and relies on manual identification of cells and
% neighbors to build connectivity matrices.
%
% Creates files for the IDs and neighbors ([image]_connections.csv) and
% coordinates ([image]_coords.csv).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mstart
load cmap
set(groot, 'defaultTextFontSize', 16)

opdir_temp = dir('output/*coords.csv');
opdir = cell(length(opdir_temp),1);
for i = 1:length(opdir_temp)
    opdir{i} = opdir_temp(i).name;
end
datdir_temp = dir('data/*.tif');
groups = {};
for i = 1:length(datdir_temp)
    thisname = strsplit(datdir_temp(i).name,'.tif');
    thisname = thisname{1};
%     if ~ismember(strjoin({thisname,'_coords.csv'},''), opdir)
        groups = vertcat(groups,thisname);
%     end
end

%%
for g = 1:length(groups)
    clear chi
    groupname = groups{g};
    disp(strjoin({num2str(g),groupname}))
    %%
    % Initialize figure
    try
        figure(fig);
    catch
        fig = figure();
    end
    clf
    fig.Name = groupname;

    % Colormap
    load cmap

    % Read in image
    [img, ~,~] = imread(strjoin({'data/',groupname,'.tif'},''));
    % [img, ~,~] = imread(strjoin({'data/',groupname,'_labels.jpg'},''));

    imshow(img)
    hold on

    %% Using existing coordinates file
    if exist(strjoin({'output/',groupname,'_coords.csv'},''),'file') == 2
        readtab = readtable(strjoin({'output/',groupname,'_coords.csv'},''));
        xax0 = readtab.x(strcmp(readtab.Label,'ax'));
        yax0 = readtab.y(strcmp(readtab.Label,'ax'));
        xap0 = readtab.x(strcmp(readtab.Label,'ap'));
        yap0 = readtab.y(strcmp(readtab.Label,'ap'));
        xtt0 = readtab.x(strcmp(readtab.Label,'tt'));
        ytt0 = readtab.y(strcmp(readtab.Label,'tt'));
        xtp0 = readtab.x(strcmp(readtab.Label,'tp'));
        ytp0 = readtab.y(strcmp(readtab.Label,'tp'));
        xbs0 = readtab.x(strcmp(readtab.Label,'bs'));
        ybs0 = readtab.y(strcmp(readtab.Label,'bs'));

        readtab = readtab(readtab.ID>0,:);
        tab = readtab(~strcmp(readtab.Label,'mc') & ~strcmp(readtab.Label,'ms'),:);
    else
        %% Select points
        fig.Name = 'Select AXIAL XYLEM (ax)';
        [xax,yax] = getpts(fig);
        plot(xax,yax,'s','Color',cmap.ax,'MarkerFaceColor',cmap.ax,'MarkerSize',15)

        fig.Name = 'Select AXIAL PHLOEM (ap)';
        [xap,yap] = getpts(fig);
        plot(xap,yap,'s','Color',cmap.ap,'MarkerFaceColor',cmap.ap,'MarkerSize',15)

        fig.Name = 'Select TRANSFUSION TRACHEIDS (tt)';
        [xtt,ytt] = getpts(fig);
        plot(xtt,ytt,'s','Color',cmap.tt,'MarkerFaceColor',cmap.tt,'MarkerSize',15)

        fig.Name = 'Select TRANSFUSION PARENCHYMA (tp)';
        [xtp,ytp] = getpts(fig);
        plot(xtp,ytp,'s','Color',cmap.tp,'MarkerFaceColor',cmap.tp,'MarkerSize',15)

        fig.Name = 'Select BUNDLE SHEATH (bs)';
        [xbs,ybs] = getpts(fig);
        plot(xbs,ybs,'s','Color',cmap.bs,'MarkerFaceColor',cmap.bs,'MarkerSize',15)

        fig.Name = groupname;

        % Combine points with IDs
        coords = [xax yax;
            xtt ytt;
            xap yap
            xtp ytp
            xbs ybs];

        coords = [coords (1:length(coords))'];

        IN = cellstr(vertcat(...
            repmat('ax',length(xax),1),...
            repmat('tt',length(xtt),1),...
            repmat('ap',length(xap),1),...
            repmat('tp',length(xtp),1),...
            repmat('bs',length(xbs),1)));

        tab = table(IN, coords(:,1), coords(:,2), coords(:,3), 'VariableNames', {'Label','x','y','ID'});

    end

    clf
    shownetwork(img,tab,cmap)

    ax = strcmp('ax',tab.Label);
    tt = strcmp('tt',tab.Label);
    ap = strcmp('ap',tab.Label);
    tp = strcmp('tp',tab.Label);
    bs = strcmp('bs',tab.Label);
    ms = strcmp('ms',tab.Label);
    air = strcmp('air',tab.Label);

    %% Edit points
    act = '';
    editpoints

    %% Select connections

    % If chi already exists, initialize with that
    if size(chi,1) ~= length(tab.ID) && exist(strjoin({'output/',groupname,'_connections.csv'},''),'file') == 2
        chi = read_connections(strjoin({'output/',groupname,'_connections.csv'},''));
    end
    editchi
    shownetwork(img,tab,cmap,chi)

    %% Make sure points are in the proper order
    ax = find(strcmp(tab.Label,'ax'));
    tt = find(strcmp(tab.Label,'tt'));
    ap = find(strcmp(tab.Label,'ap'));
    tp = find(strcmp(tab.Label,'tp'));
    bs = find(strcmp(tab.Label,'bs'));
    %     ms = find(strcmp(tab.Label,'ms'));
    %     air = find(strcmp(tab.Label,'air'));
    %     bx = find(strcmp(tab.Label,'bx'));
    %     bp = find(strcmp(tab.Label,'bp'));

    order = [ax;tt;ap;tp;bs];
    tab = tab(order,:);
    N = sum(ismember(tab.Label, {'ax','tt','ap','tp','bs'}));
    tab.ID(1:N) = (1:N)';
    chi = chi(order,order);


    %% Add boundary xylem & phloem

    if sum(strcmp(tab.Label,'bp')) == 0
        fig.Name = 'Select BOUNDARY PHLOEM (bp)';
        [xbp,ybp] = getpts(fig);
        tab((end+1),:) = table({'bp'},xbp(1),ybp(1),0);
        plot(xbp,ybp,'s','Color',cmap.bp,'MarkerFaceColor',cmap.bp,'MarkerSize',15)
    end

    if sum(strcmp(tab.Label,'bx')) == 0
        fig.Name = 'Select BOUNDARY XYLEM (bx)';
        [xbx,ybx] = getpts(fig);
        tab((end+1),:) = table({'bx'},xbx(1),ybx(1),-1);
        plot(xbx,ybx,'s','Color',cmap.bx,'MarkerFaceColor',cmap.bx,'MarkerSize',15)
    end

    fig.Name = groupname;
    figure(fig)
    shownetwork(img,tab,cmap,chi);


    %% Add air & mesophyll

    theta = pi/8;

    ax = strcmp('ax',tab.Label);
    tt = strcmp('tt',tab.Label);
    ap = strcmp('ap',tab.Label);
    tp = strcmp('tp',tab.Label);
    bs = strcmp('bs',tab.Label);
    ms = strcmp('ms',tab.Label);
    air = strcmp('air',tab.Label);
    cpairs = get_pairs(chi,tab);
    if sum(ms)~=sum(bs)
        if sum(air)>0
            tab(strcmp('air',tab.Label),:) = [];
        end
        if sum(ms)>0
            tab(strcmp('ms',tab.Label),:) = [];
        end
        L = length(tab.Label);
        N = sum(ismember(tab.Label,{'ax','tt','ap','tp','bs'}));
        nbs = sum(bs);
        ms_str = mat2cell(repelem('ms',nbs,1),nbs,2);

        BS = find(bs);

        ravg = mean((range(tab.x(cpairs),2).^2 + range(tab.y(cpairs),2).^2).^0.5);

        rcent = ([tab.x(BS(1)) tab.y(BS(1))] + [tab.x(BS(end)) tab.y(BS(end))])/2;

        for i = 1:nbs
            r = [tab.x(BS(i)) tab.y(BS(i))] - rcent;
            r = ravg*r/norm(r);
            mscoords_i = ([tab.x(BS(i)) tab.y(BS(i))]' + ...
                [cos(theta) -sin(theta); sin(theta) cos(theta)]*r')';
            aircoords_i = ([tab.x(BS(i)) tab.y(BS(i))]' + ...
                2*[cos(theta/3) sin(theta/3); -sin(theta/3) cos(theta/3)]*r')';
            newcoords_i = [mscoords_i; aircoords_i];
            tab(L+[0 nbs]+i,:) = table({'ms';'air'},newcoords_i(:,1),newcoords_i(:,2),[N+i; -200-N+nbs-i]);
        end
    end

    shownetwork(img,tab,cmap,chi)

    editpoints


    %% write coords and connections to csv file and an image of the network
    % Coords
    writetable(tab, strjoin({'output/',groupname,'_coords.csv'},''))

    % Connections
    ncol = max(sum(chi,2));
    fid = fopen(strjoin({'output/',groupname,'_connections.csv'},''),'w');

    for i = 1:length(tab.Label)
        if ismember(tab.Label{i}, {'bx','bp','air','ms'})
            continue
        end
        fprintf(fid,'%s',tab.Label{i});
        fprintf(fid,',%d',find(chi(i,:)==1));
        fprintf(fid,repelem(',',ncol-sum(chi(i,:))));
        if i < length(tab.Label)
            fprintf(fid,'\n');
        end
    end
    fclose(fid);

    % Save image to images/
    saveas(fig,strjoin({'images/',groupname,'.png'},''))

end