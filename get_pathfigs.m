%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function get_pathfigs
% Melissa H Mai
%
% Generate a scatterplot of each cell in the network as a function of its
% pathlength from a reference node.
%
% INPUTS
% config        Configuration name, as the prefix for all the relevant
%                   output file [config].mat. NOTE: tt_model must be run
%                   with 'writemat', true.
%
% xvar          Name of the x-variable, as given in the results .mat file
% 
% yvar          Name of the y-variable, as given in the results .mat file
%
% showlegend    Boolean value whether to show the legend
%
%
% OUTPUT
% thisax        Axis handle for the resulting figure
%
% leg           Handle for the resulting legend
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [thisax,leg] = get_pathfigs(config,xvar,yvar,showlegend)
    if ~exist('showlegend','var'); showlegend = false; end
    mstart
    load('cmap.mat');
    cells = {'tt','tp','ax','ap','en','ms'};
    shapes = struct('ax','^','ap','v','tt','d','tp','o','en','s','ms','d');
    env = load(strjoin({'output/',config,'.mat'},''));
    fnames = fieldnames(env);
    run = env.(fnames{1});
    
    data = needle_analysis(run,true);
    data.ax = data.ex;
    data.ap = data.ep;
    
    x = data.(xvar);
    y = data.(yvar);
    
    thisax = gca;
    for c = 1:length(cells)
        if (length(x) > length(y)) && (c == length(cells))
            continue
        end
        p(c)=plot(x(data.(cells{c})), y(data.(cells{c})),shapes.(cells{c}),...
            'MarkerFaceColor',cmap.(cells{c}), 'MarkerEdgeColor',cmap.(cells{c})*.7,...
            'Color',cmap.(cells{c}),'DisplayName',cells{c},'MarkerSize',8);
        hold on
    end
    
    xlim([min(x)-0.5 max(x)+0.5])
    set(gca,'FontSize',14)
    xlabel(xvar)
    ylabel(yvar)

    if showlegend
        leg = legend(p([3 4 1 2 5:end]), cells([3 4 1 2 5:end]));
    else
        leg = [];
    end
    
    hold off
end