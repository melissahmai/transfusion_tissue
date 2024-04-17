%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function shownetwork
% Melissa H Mai
%
% Display the current network during network construction (build_network)
%
% INPUT
% img       Image handle
% tab       Coordinate table with cell IDs
% cmap      Colormap (from cmap.mat)
% chi       Connectivity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function shownetwork(img,tab,cmap,chi)
    ax = strcmp('ax',tab.Label);
    tt = strcmp('tt',tab.Label);
    ap = strcmp('ap',tab.Label);
    tp = strcmp('tp',tab.Label);
    bs = strcmp('bs',tab.Label);
    ms = strcmp('ms',tab.Label);
    air = strcmp('air',tab.Label);
    bp = strcmp('bp',tab.Label);
    bx = strcmp('bx',tab.Label);
    nbs = sum(bs);

    hold off
    imshow(img)
    hold on
    if exist('chi','var') 
        for i = 1:(length(chi))
            cnx = find(chi(i,:)==1);
            if isempty(cnx)
                continue
            end
            for j = 1:length(cnx)
                if any(ismember(tab.Label(ismember(tab.ID,[i cnx(j)])), {'ax','tt'}))
                    color = 'blue';
                else
                    color = 'red';
                end
                pts = find(ismember(tab.ID,[i cnx(j)]));
                plot(tab.x(pts), tab.y(pts), '-', 'Color', color)
            end
        end
    end
    
    if sum(ms)>0
        ms_id = find(ms);
        bs_id = find(bs);
        for i = 1:length(ms_id)
            if sum(air)>0
                air_id = find(tab.ID==(-200-tab.ID(bs_id(i))));
                pathx = tab.x([air_id ms_id(i) bs_id(i)]);
                pathy = tab.y([air_id ms_id(i) bs_id(i)]);
            else
                pathx = tab.x([ms_id(i) bs_id(i)]);
                pathy = tab.y([ms_id(i) bs_id(i)]);
            end
            plot(pathx,pathy, ...
                '-', 'Color', 'green') 
        end
    end
    
    if sum(bp)>0
        plot(tab.x(logical(bp+ap)),tab.y(logical(bp+ap)),'-','Color','green')
    end
    if sum(bx)>0
        plot(tab.x(logical(bx+ax)),tab.y(logical(bx+ax)),'-','Color','green')
    end
    
    plot(tab.x(ax),tab.y(ax),'s','Color',cmap.ax,'MarkerFaceColor',cmap.ax,'MarkerSize',15)
    plot(tab.x(tt),tab.y(tt),'s','Color',cmap.tt,'MarkerFaceColor',cmap.tt,'MarkerSize',15)
    plot(tab.x(ap),tab.y(ap),'s','Color',cmap.ap,'MarkerFaceColor',cmap.ap,'MarkerSize',15)
    plot(tab.x(tp),tab.y(tp),'s','Color',cmap.tp,'MarkerFaceColor',cmap.tp,'MarkerSize',15)
    plot(tab.x(bs),tab.y(bs),'s','Color',cmap.bs,'MarkerFaceColor',cmap.bs,'MarkerSize',15)
    plot(tab.x(ms),tab.y(ms),'s','Color',cmap.ms,'MarkerFaceColor',cmap.ms,'MarkerSize',15)
    plot(tab.x(air),tab.y(air),'s','Color',cmap.air,'MarkerFaceColor',cmap.air,'MarkerSize',15)
    plot(tab.x(bp),tab.y(bp),'s','Color',cmap.bp,'MarkerFaceColor',cmap.bp,'MarkerSize',15)
    plot(tab.x(bx),tab.y(bx),'s','Color',cmap.bx,'MarkerFaceColor',cmap.bx,'MarkerSize',15)
    text(tab.x, tab.y, num2str(tab.ID),'HorizontalAlignment','Center')
    
    
    xrange = range(tab.x);
    yrange = range(tab.y);
    xlim([min(tab.x)-xrange/10 max(tab.x)+xrange/10])
    ylim([min(tab.y)-yrange/10 max(tab.y)+yrange/10])

end