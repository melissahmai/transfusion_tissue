%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function needle_analysis
% Melissa H Mai
%
% Takes the model output from tt_model and extracts relevant information,
% and can also convert values to SI units.
%
% INPUTS
% modelOp           Output from tt_model
% convert           Boolean whether to convert to SI units
%
% OUTPUT
% Structure     with the following fields:
% note          Note inherited from tt_model
% k             Number of reported steps
% Time          Final time to end of simulation
% C             Vector of concentrations for each node
% S             Vector of starch for each node
% Cp            Final concentration in the phloem
% P             Vector of pressures for each node
% Psi           Vector of water potentials for each node
% Cm            Vector of mesophyll concentrations
% RT            Thermal energy
% Ux            Water flow from the xylem (>0 for proper function)
% Up            Water flow into the phloem (<0 for proper function)
% E             Total evaporation (<0 for proper function)
% J             Sugar export through the phloem (>0 for proper function)
% wue           Export efficiency Jp/E
% recycling     Fraction of water from the xylem returned to the phloem
% sugar_in      Sugar flux into the BS from MS
% sugar_back    Sugar backwashing into the mesophyll
% backflow      Ratio of sugar flow backwashing into the mesophyll   
% ax, ex, ap,   Boolean vectors that indicate cell identities
%   ep, tt,
%   tp, bs,
%   ms
% paths_xy      Shortest pathlengths to each node from the axial xylem
% paths_ph      Shortest pathlengths to each node from the axial phloem
% epath_xy      Shortest pathlengths to each connection from the axial xylem
% epath_ph      Shortest pathlengths to each connection from the axial phloem
% epath_ms      Shortest pathlengths to each connection from the axial mesophyll
% Pe            Matrix of Peclet numbers for each connection
% coflow        For each discrete pathlength from the xylem, the proportion of
%                   connections in which water and sugar travel in the same direction
% trachflow     For each cell, water flow contributions from tracheary elements
% Ui            Positive component of water fluxes through each node
% Ji            Positive component of sugar fluxes through each node
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function op = needle_analysis(modelOp, convert)
    % Unit conversion, if specified
    if exist('convert', 'var') && convert
        convert = true;
        modelOp.C = [modelOp.C(end,:)'; modelOp.Cm(end,:)'];
        modelOp = modelOp.SIconvert(modelOp);
    else
        convert = false;
    end
    
    % Unpack information from passed structure
    unpackStruct(modelOp)

    op = struct;

    % Carry over notes
    op.note = note;
    
    %% Take last time point
    op.k = k;
    if ismember('args',fieldnames(modelOp)) && ismember('k0',fieldnames(modelOp.args))
        op.Time = Time(length(Time)-modelOp.args.k0+1);
    else
        op.Time = Time(length(Time)-49);
    end
    if ~convert
        op.C = [C(end,:)'; Cm(end,:)'];
    else 
        op.C = C;
    end
    op.S = S(end,:);
    op.Sn = op.S/mean(op.S(op.S>0));
    op.Cp = Cp(end);
    op.P = [P(end,:)'; Pm(end,:)'];
    op.Psi = op.P - RT*op.C;
    op.Cm = Cm(end,:)';
    op.RT = RT;

    %% Rescaling by Lp
    Uib = Uib*lpscale;
    Jib = Jib*lpscale;
    Jij = Jij*lpscale;
    Uij = Uij*lpscale;
    
    %% Basic/overview metrics
    % Flow into the xylem
    op.Ux = sum(Uib(strcmp(IN,'ax') | strcmp(IN,'ex')));

    % Flow out of the phloem
    op.Up = sum(Uib(strcmp(IN,'ap') | strcmp(IN,'ep')));

    % Evaporative flux
    op.E = -(op.Ux+op.Up);

    % Sugar Export
    op.J = -sum(Jib(strcmp(IN,'ap') | strcmp(IN,'ep')));

    % Water Use Efficiency metrics
    op.wue = -op.J/op.E;
    op.recycling = -op.Up/op.Ux;

    % Backflow of sugar from the bundlesheath
    % Get sugar flux through the mesophyll cells
    Jib_meso = sum(Jij((length(IN)+1):end,:),2);
    
    % Negative flux means sugar is leaving the mesophyll into BS
    op.sugar_in = sum(-Jib_meso(Jib_meso < 0));
    
    % Positive flux means sugar is entering the mesophyll (backflow)
    op.sugar_back = sum(Jib_meso(Jib_meso > 0));
    
    % Get the ratio of sugar going back over sugar going in
    op.backflow = op.sugar_back/op.sugar_in;
    
    %% IDs
    N = length(IN);
    ids = {'ax','ex','ap','ep','tt','tp','bs','ms'};
    for i = 1:length(ids)
        op.(ids{i}) = strcmp(IN, ids{i});
    end
    op.ax = op.ex | strcmp(IN,'ax');
    op.ap = op.ep | strcmp(IN,'ap');
    op.ms = (length(IN)+1):size(chi,1);
    
    %% Paths

    % To exclude tt in path finding
    chi_nott = chi;
    chi_nott(op.tt,:) = 0;
    chi_nott(:,op.tt) = 0;
    
    % But connect the tp that's attached to ap to the ax
    chi_nott(op.ex | op.ax,chi_nott((op.ep | op.ap),:)==1) = 1;
    chi_nott = chi_nott - diag(chi_nott).*eye(size(chi_nott,1));
    op.paths_xy = findpath(chi_nott,1);
    op.paths_ph = findpath(chi_nott,(op.ep | op.ap));
    op.epath_xy = findpath_edge(chi_nott,1);
    op.epath_ph = findpath_edge(chi_nott,(op.ep | op.ap));
    op.epath_ms = findpath_edge(chi_nott,op.ms);
    
    %% Detail metrics

    % Pe of each connection
    op.Pe = Uc./Dc;
    
    % U and J in the same direction
    coflow = sign(Uij.*Jij);
    coflowpath = op.epath_ms(coflow ~= 0);
    coflow = coflow(coflow~=0);
    [coflowpath,I] = sort(coflowpath);
    coflow = coflow(I);
    paths = unique(coflowpath);
    op.coflow = zeros(length(paths),4);
    for i = 1:length(paths)
        npath = sum(coflowpath == paths(i));
        nco = sum(coflow(coflowpath==paths(i))==1);
        op.coflow(i,:) = [paths(i), nco, npath-nco, nco/npath];
    end
    
    % Tracheid flow
    op.trachflow = sum(Uij(:,strcmp(IN,'ax')|strcmp(IN,'ex')|strcmp(IN,'tt')),2);
    op.trachflow(strcmp(IN,'ax')|strcmp(IN,'ex')) = ...
        op.trachflow(strcmp(IN,'ax')|strcmp(IN,'ex')) +...
        Uib(strcmp(IN,'ax')|strcmp(IN,'ex'));
    
    % Total flow
    Up = Uij.*(Uij>0);
    op.Ui = sum(Up,2) + Uib.*(Uib>0);

    Jp = -sum(Jij.*(Jij<0),2);
    Jp(1:N) = Jp(1:N) - Jib(1:N).*(Jib(1:N)<0);
    op.Ji = Jp;
end