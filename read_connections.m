%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function read_connections
% Melissa H Mai
%
% Reads in connection information from [config]_connections.csv to prepare
% inputs for tt_model
%
% INPUTS
% filename      Filename (eg [config]_connections.csv)
% axap          Boolean whether to override and enforce an ax-ap connection
%                   Default: false
%
% OUTPUTS: chi, IN      See tt_model for description
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [chi,IN] = read_connections(filename, axap)
    if ~exist('axap','var'); axap = false; end
    
    % Read in file
    opts = detectImportOptions(filename);
    ncol = length(opts.VariableNames);
    format = strjoin({'%s',repmat('%f',1,ncol-1)},'');
    fid = fopen(filename);
    connections = textscan(fid,format,'Delimiter',',');
    fclose(fid);
    
    % Get node identifier strings
    IN = connections{1};
    
    if isempty(IN{end})
        N = length(IN)-1;
    else
        N = length(IN);
    end
    IN = IN(1:N);
    
    % Connectivity matrix
    connections = cell2mat(connections(2:end));
    connections = connections(1:N,:);
    

    chi = zeros(N);
    for i = 1:(N-1)
        chi(i,connections(i,~isnan(connections(i,:)))) = 1;
    end
    
    if axap
        chi(1, strcmp(IN,'ap')) = 1;
    end
end