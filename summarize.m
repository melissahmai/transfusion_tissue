%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function summarize
% Melissa H Mai
%
% Summarize model analyses (from needle_analysis)
%
% INPUT
% Tables of node identities and coordinates (from [config]_coords.csv). Multiple 
%       tables can be passed as individual arguments.
%
% OUTPUT
% tab       Summary table of the following fields:
%               Run: run name
%               Cm: Average mesophyll concentration
%               Cp: Phloem concentration
%               Pph: Phloem pressure
%               Pm: Mesophyll pressure
%               Ux: Xylem water flow
%               Up: Phloem water flow
%               E: Transpiration
%               J: Sugar export
%               WUE: Export efficiency
%               Recycling: Phloem water flow / Xylem water flow
%               Sugar_in: Sugar going into BS from MS
%               Sugar_back: Sugar backwashing into MS
%               Backflow: Sugar_back/Sugar_in
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tab = summarize(varargin)
    tab = table;
    k = 1;
    for i = 1:length(varargin)
        this = varargin{i};
        for j = 1:length(this)
            now = this{j};
            P_ph = now.P(now.ep);
            P_m = mean(now.P(now.ms));
            tab(k,:) = cell2table({inputname(i), mean(now.Cm), now.Cp, P_ph, P_m, now.Ux, now.Up, now.E,...
                now.J, now.wue, now.recycling, ...
                now.sugar_in, now.sugar_back, now.backflow});
            k = k+1;
        end
    end
    
    tab.Properties.VariableNames = {'Run', 'Cm','Cp', 'Pph', 'Pm', 'Ux', 'Up', 'E', 'J', 'WUE', ...
        'Recycling', 'Sugar_in', 'Sugar_back', 'Backflow'};
end