%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function get_configs
% Melissa H Mai
%
% List all available configurations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function configs = get_configs
    Chi=struct;IN=struct;%runs=struct;
    
    % First pass: default values (ie, constant mesophyll concentration)
    % Second pass: constant mesophyll flux boundary condition
    configdir_temp = dir('output/*connections.csv');
    configs = cell(length(configdir_temp),1);
    for i = 1:length(configs)
        thisname = strsplit(configdir_temp(i).name, '_connections');
        configs{i} = thisname{1};
    end
end