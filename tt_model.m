%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function tt_model
% Melissa H Mai
%
% Main function for the transfusion tissue network model.
% Related scripts/functions:
%       build_network: construct networks from images
%       read_connections: get IN and chi inputs
%       needle_analysis: perform basic analysis of model output
%       summarize: compile analyses of multiple runs
%
% INPUTS (mandatory)
% IN        Identity names of the N-many nodes, in order. Identity codes are
%                   'ex' : edge axial xylem (directly connected to boundary)
%                   'ax' : axial xylem
%                   'tt' : transfusion tracheid
%                   'ep' : edge axial phloem (directly connected to boundary)
%                   'ap' : axial phloem
%                   'tp' : transfusion parenchyma
%                   'bs' : bundle sheath cell
%               All nodes of a type must be grouped together, in the order
%               of 'ax'/'ex', 'tt', 'ap'/'ep', 'tp', 'bs'.
%               If any node is labeled as 'ex', then only that node will be
%               treated as a boundary-adjacent node with the xylem condition.
%               All other 'ax' will be treated as completely internal nodes (ie,
%               similar to the 'tt' nodes). If no 'ex' is specified, any 'ax'
%               will then be treated as boundary nodes. The same goes for 'ep'
%               and the phloem nodes.
%
% chi       The connection matrix chi is a NxN boolean matrix, preferably given
%               as an upper triangular matrix (to avoid errors in symmetry
%               from manual definition of connections). The identity line
%               should be zero. Each cell at position (i,j) corresponds to
%               a potential connection between cells i and j. If there is a
%               connection, chi(i,j) = 1, otherwise it is zero. Only the
%               upper triangle of the inputted matrix will be read.
%
%
% INPUTS (optional)
% Optional inputs should be passed as 'Name', Value pairs when calling the
% function (ie, tt_model(IN, chi, 'c0', 0, 'nsteps', 1000). The first optional
% input may be passed as either a structure (as in the .args output from a 
% previous simulation) or a cell array, followed by additional Name/Value pairs.
% ie,   tt_model(IN, chi, struct('c0', 0, 'nsteps', 1000), 'tau', 1e-2)
%   or  tt_model(IN, chi, prev_run.args, 'tau', 1e-2)
%   or  tt_model(IN, chi, {'c0', 0, 'nsteps', 1000}, 'tau', 1e-2)
%
% 'nsteps', 'nstep'     Maximum number of steps
%                           Default: 50000
%
% 'rdt'                 Readout increment. Only every rdt-th step will be
%                           recorded to reduce memory strain.
%                           Default: 10
%
% 'k0', 'kmin', 'mink'  Minimum number of recorded steps. This is specified
%                           to help ensure steady-state is established.
%                           Default: 50
%
% 'c0', 'cinit'         Nx1 vector of initial conditions for concentrations,
%                           where c0(i) corresponds to node i, with identity and
%                           connections defined in IN and chi, respectively.
%                           It is possible for c0 to be an N_lx1 matrix, where
%                           N_1 is the number of living cell nodes (ie, 
%                           excluding the N_t tracheary elements). In this case,
%                           c0(i) corresponds to node N_t+i.
%                           Default: All concentrations 0
%
% 'cms', 'c_ms',        Concentration of sugar in the mesophyll (sugar
%   'cm', 'c_m'             producing boundary) cell, passed either as a single
%                           value or as an N_m x 1 vector. This serves as the
%                           fundamental concentration for dimensional
%                           analysis. If constant sugar flux is used instead, cm
%                           may be set to 0; in this case, the value for
%                           dimensional analysis will still be the default
%                           value. Otherwise, the value used will be the mean of
%                           the mesophyll cells.
%                           Default: 500
%
% 'cm_star', 'cmstar'   Reference mesophyll concentration, if the mesophyll
%                           concentation is not uniform or constant.
%                           Default: 500 (if cm not specified), or
%                                    mean(cm)
%
% 'Jm', 'Jms', 'J_m',   Fixed sugar production rate, per mesophyll cell.
%   'J_ms'                  At each step, Jm*dt sugar is added to each
%                           mesophyll cell. This *does not* specify the
%                           sugar flux from each mesophyll cell to its
%                           corresponding bundle sheath cell, but rather
%                           the photosynthetic rate. When Jm = 0, sugar is
%                           not produced, but any residual concentration of
%                           sugar is drained from the mesophyll.
%                           Default: unspecified (cm is constant)
%
% 'allow_buildup',      Boolean value allowing the concentration of the 
%   'buildup'               mesophyll to exceed its set value cm if sugar
%                           is allowed to backwash from the bundle sheath.
%                           Default: true
%
% 'maxcm', 'cm_max',    Maximum factor cm is allowed to exceed its set
%   'cmmax',                value cm, if allow_buildup is true.
%   'buildup_max'           Default: 1.05
%
% 'c_b', 'c_bounds',    Concentrations at the boundary conditions, in the
%   'cb', 'cbounds'         order and scaled to cms
%                           [xylem phloem], or [phloem]
%
%                           If only one quantity is provided, it is assumed
%                           that the xylem has a concentration of 0.
%
%                           If the concentration given for the phloem exceeds 1,
%                           it is assumed that the phloem has not been scaled to
%                           the mesophyll concentration and will then be
%                           rescaled by the mesophyll concentration.
%                           Default: [0; 0]
%
% 'E', 'ET'             Transpiration rate. If unspecified, the boundary
%                           condition at the bundle sheath and mesophyll
%                           will be a constant air vapor pressure, with
%                           variable transpiration.
%
% 'EperBS', 'perBS'     Boolean defining whether the specified
%                           transpiration E is applied over the whole
%                           needle or scaled by the number of bundle sheath
%                           cells.
%                           Default: true
%
% 'P_b', 'P_bounds',    Pressure at the boundary conditions, in the order
%   'Pb', 'Pbounds'         [xylem phloem air]
%                           to be given in units of bars (ie 1e5 Pa = 1).
%                           P_air is overridden if the transpiration E is
%                           specified.
%                           Default: [0;1;-15]
%
% 'npd', 'n_pd'         Permeabilities at connection interfaces
%   'lp', 'lpstruct'        If given as a struct, it can determine either
%                           group or individual densities. To define entire
%                           groups, use the fields
%                               'trach' for tracheary connections
%                               'osm'   for osmotic connections
%                               'pd'    for plasmodesmal connections
%
%                           Individual connections can also be defined
%                           using the fields (and corresponding to the
%                           groups denoted in parentheses)
%                               (trach) : 'axax','axtt','tttt'
%                               (osm)   : 'axap','axtp','ttap','tttp',
%                                         'axbs','ttbs'
%                               (pd)    : 'apap','aptp','tptp','apbs',
%                                         'tpbs','bsbs','bsms'
%
%                           A combination of group and individual
%                           definitions can be used, in which case
%                           individual definitions are given priority and
%                           unspecified connections are defined by the group.
%                           If neither the individual nor group value is
%                           given for a connection, it will be assigned a
%                           default value.
%
%                           Otherwise, npd can be passed as a 16x1 vector,
%                           in the order 
%                               'axax','axtt','tttt','axap','axtp','ttap',
%                               'tttp','axbs','ttbs','apap','aptp','tptp',
%                               'apbs','tpbs','bsbs','bsms'
%
%                           If only a single value is passed, all Lp will
%                           be set equal to each other.
%
%                           The final values will then be scaled so that
%                           the largest value is 1. This is ok since all
%                           terms are multiplied by some form of the
%                           permeability except for the sugar flux j, in
%                           which case this scale is absorbed into the time
%                           step tau. As a result, npd (number of plasmodesmata)
%                           is considered equivalent to the permeability in this
%                           situation.
%                        
%                           Defaults: (trach) : 1
%                                     (osm)   : 0.1
%                                     (pd)    : 0.1
%
% 'scalelp',            Boolean value stating whether the permeabilities as 
%   'scale_lp',             given in lp (see above) should be rescaled so the
%   'lpscale', 'lp_scale'   highest value of LP is 1. NOTE: if only certain
%                           connections are specified in lp (above) and if
%                           scalelp is true, then ALL LP values will be rescaled
%                           (including the default ones). Therefore, it is
%                           advised that scalelp = false if only certain LP
%                           values are specified
%                           Default: true
%                  
% 'lp_bounds', 'lp_b',  Permeabilities characterizing the connections to
%   'lp_bound', 'lpb',      boundaries, given in the order
%   'lpbounds',             [xylem phloem air]
%   'lpbound'               Default: [LP_axtt*.02; LP_aptp*.2; max(LP)*1.5e-3]'
%
% 'bsair', 'airbs',     Boolean value stating whether water flow is impeded
%   'bsa', 'casparian'      between the bundle sheath cells and the air spaces
%                           Default: false (ie, flow allowed between bs & air)
%
% 'sig', 'sigma'        Reflection coefficient, 0 <= sig <= 1. 
%                           Default: calculated using the equation
%                           found in Deen & Dechadilok 2006,
%                           approx 0.3
%
% 'psig', 'phloemsigma' Boolean value defining whether the connection
%   'psigma'                between the loading and boundary phloem is
%                           characterized by sig = 0 (false) or sig as
%                           defined above (true)
%                           Default: false
%
% 'omega', 'w'          Effective diffusion omega through the PD. If no value is
%                           given, omega will be set according to the
%                           calculation for the hindrance factor H0 (see code).
%                           NOTE: omega will be scaled according to scalelp (see
%                           above), in the same way the Lp will be, if scalelp
%                           is true.
%                           Default: 3*H0*Lp
%
% 'rs', 'rsolute',      Hydrodynamic radius of the solute
%   'r_s', 'r_solute'       Default: 4.2e-10
%
% 'T', 'temp',          Temperature, in Kelvin
%   'temperature'           Default: 298.15
%
% 'dt', 'tau'           Initial relaxation/update step size used for sugar
%                           calculation. Will be replaced during adaptive time
%                           stepping
%                           Default: 5e-2
%
% 'epsilon', 'eps',     Tolerance for convergence criterion such that when
%   'tol', 'tolerance'      norm(all concentration changes) <= eps and the
%                           change in phloem concentration over k0 steps (see
%                           below) divided by the current time step tau <= eps,
%                           the simulation will stop. Otherwise, the
%                           simulation will end when nsteps has been
%                           reached to avoid infinite computation. If,
%                           instead, you wish to run the simulation to
%                           nsteps, set epsilon = 0. This epsilon is also used
%                           for updating the adaptive time step; in particular,
%                           twice of this epsilon value is used to determine the
%                           accuracy of a given time step.
%                           Default: 1e-6
%
% 'h'                   Width of the plasmodesmal sheath; this serves as the 
%                           characteristic length for dimensional analysis.
%                           Default: 1e-9
%
% 'D', 'Dcyt', 'Dfree'  Cytosolic/free diffusion coefficient for the solute
%                           (sugar). This serves as a primary variable for
%                           dimensional analysis.
%                           Default: 2.3e-10
%
% 'eta', 'visc',        Viscosity of the medium. This serves as a primary
%   'viscosity'             variable for dimensional analysis
%                           Default: 2e-3
%
% 'ks', 'kstarch'       Rate constants for starch kinetics. Given as a structure
%   'k_s', 'k_starch'       with the fields 
%                               'vmax': Michaelis-Menten Vmax
%                               'km'  : Michaelis constant
%                               'kh'  : Hydrolysis/digestion rate constant.
%                           If ks is passed as a single value, all rates will be
%                           set to that value (eg, 0).
%                           Default: struct('vmax',1,'km',1,'kh',1)
%
% 's0', 'sinit'         Nx1 vector of initial conditions for starch 
%                           concentrations, where c0(i) corresponds to node i.
%                           See c0 above.
%                           Default: All 0
%
% 'filename',           Filename of the desired output tables. The tables
%   'writeto'               will be written to [filename]_nodes.csv and
%                           [filename]_edges.csv. All input arguments will be
%                           saved as [filename].mat, which will include
%                           all outputs from the model. All files will be
%                           written to the output/ directory.
%                           If not given, no file will be written.
%
% 'writemat',           Boolean defining whether the output structure should be
%   'savemat'               saved to a file [filename].mat in a structure called
%                           [filename]_run. Saving these structures can take up
%                           a lot of storage space.
%                           Default: false
%
% 'note', 'notes'       Any notes to include with the output structure for
%                           future reference.
%
%
% OUTPUTS
% The output is a structure `output` with the following fields:
% C         Matrix of the concentrations at each node, with each row
%               corresponding to one time step. The initial condition is
%               not included. Given in terms of the simulation units.
%
% P         Matrix of pressures at each node. Same properties as C above.
%
% S         Matrix of starch content at each node. Same properties as C above.
%
% Psi       Water potential at each node. Same properties as C and P above.
%
% Uij       Anti-symmetric matrix of water fluxes between nodes at steady 
%               state. The identity line is zero.
%
% Jij       Anti-symmetric matrix of sugar fluxes between nodes. See Uij.
%
% Uc        Anti-symmetric matrix of advective sugar fluxes between nodes. See
%               Uij and Jij.
%
% Dc        Anti-symmetric matrix of diffusive sugar fluxes between nodes. See
%               above. Uc + Dc = Jij
%
% Uib       Water flows at the boundaries at the final solution
%
% Jib       Sugar flows at the boundaries at the final solution
%
% k         Number of steps used to finish simulation
%
% Time      Vector of time at each readout step
%
% Cp        Evolution of the phloem concentration. Does not include the initial
%               condition.
%
% Pm        Evolution of the mesophyll pressure
%
% Cm        Evolution of the mesophyll sugar concentration
%
% Jp        Evolution of the sugar flux through the phloem
%
% ks        Rate constants for starch kinetics
%
% RT        RT in simulation units
% 
% sigma     Reflection coefficient
%
% IN        String vector of node identities
%
% chi       Boolean connectivity matrix
%
% Chi       Connectivity matrix with connection type denoted by a unique code:
%               0: none     4: ax-ap    8: ax-bs    12: tp-tp
%               1: ax-ax    5: ax-tp    9: tt-bs    13: ap-bs
%               2: ax-tt    6: tt-ap    10: ap-ap   14: tp-bs
%               3: tt-tt    7: tt-tp    11: ap-tp   15: bs-bs
%               
% Lp        16x1 vector of permeabilities/conductances. See lp in inputs above.
%               Values will be scaled if scalelp is true
%
% Lp_bounds 3x1 vector of permeabilities of boundary connections. See lp_bounds
%               above. Values will be scaled if scalelp is true
%
% lpscale   If scalelp is true, lpscale will give the scale factor used for
%               scaling the permeabilities (1 if scalelp is false, or if already
%               properly scaled). To recover original values, multiply Lp by
%               lpscale
%
% args      All optional arguments that were originally passed to tt_model
%
% note      Whatever note was passed to the original function call
%
% SIconvert Function used to convert values from this output structure to SI
%               units. NOTE: not all units will be converted (like Lp). See
%               toSI.m for more information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = tt_model(IN, chi, varargin)
    %% Input Processing
    args = struct();
    
    % If first optional input is passed as a structure or cell of arguments,
    % expand and add to whatever other arguments are passed afterwards
    if ~isempty(varargin)
        if isstruct(varargin{1})
            temp = {};
            fn = fieldnames(varargin{1});
            for i = 1:length(fn)
                temp((2*i-1):(2*i)) = {fn{i}, varargin{1}.(fn{i})};
            end
            varargin = horzcat(temp, varargin(2:end));
        elseif iscell(varargin{1})
            temp = varargin{1};
            varargin = horzcat(temp, varargin(2:end));
        end
        
        for l = 1:2:length(varargin)
            this = varargin{l+1};
            
            % Store in a structure for reporting in the output
            args.(varargin{l}) = this;
            switch lower(varargin{l})
                case {'c0', 'cinit'}
                    c0 = this;
                case {'h'}
                    h = this;
                case {'d', 'dcyt', 'dfree'}
                    D = this;
                case {'eta', 'visc', 'viscosity'}
                    eta = this;
                case {'cm', 'cms', 'c_m', 'c_ms'}
                    cm = this;
                case {'cmstar', 'cm_star'}
                    cm_star = this;
                case {'jm', 'jms', 'j_m', 'j_ms'}
                    jm = this;
                case {'rs', 'r_s', 'rsolute', 'r_solute'}
                    rs = this;
                case {'t', 'temp', 'temperature'}
                    T = this;
                case {'dt', 'tau'}
                    tau = this;
                case {'nsteps', 'nstep'}
                    nsteps = this;
                case {'epsilon', 'eps', 'tol', 'tolerance'}
                    epsilon = this;
                case {'pb', 'p_b', 'pbounds', 'p_bounds'}
                    Pb = this;
                case {'cb', 'c_b', 'cbounds', 'c_bounds'}
                    cb = this;
                case {'allow_buildup', 'buildup'}
                    allow_buildup = this;
                case {'cmmax', 'maxcm', 'cm_max','buildup_max'}
                    maxcm = this;
                case {'lp_bounds', 'lp_bound', 'lpbounds', 'lpbound', 'lpb', 'lp_b'}
                    Lp_bounds = this;
                case {'sigma', 'sig'}
                    sig = this;
                case {'lp', 'lpstruct', 'npd', 'n_pd'}
                    LP = this;
                case {'e','et'}
                    E = this;
                case {'eperbs', 'perbs'}
                    EperBS = this;
                case {'phloemsigma', 'psigma', 'psig'}
                    psig = this;
                case {'filename', 'writeto'}
                    filename = this;
                case {'writemat', 'savemat'}
                    writemat = this;
                case {'casparian', 'bsair', 'airbs', 'bsa'}
                    bsair = this;
                case {'scalelp', 'scale_lp', 'lpscale', 'lp_scale'}
                    scale_LP = this;
                case {'omega', 'w'}
                    omega = this;
                case {'ks', 'k_s', 'kstarch', 'k_starch'}
                    ks = this;
                case {'s0', 'sinit'}
                    s0 = this;
                case {'k0','kmin','mink'}
                    k0 = this;
                case {'note', 'notes'}
                    note = this;
                case {'rdt'}
                    rdt = this;
                otherwise
                    warning(strjoin({'Input', varargin{l}, 'not recognized.'}))
            end
        end
    end
    
    % Default values
    if ~exist('h','var'); h = 1e-9; end
    if ~exist('D','var'); D = 2.3e-10; end
    if ~exist('eta','var'); eta = 2e-3; end
    if ~exist('cm','var'); cm = 500; end
    if ~exist('allow_buildup','var'); allow_buildup = true; end
    if ~exist('maxcm','var'); maxcm = 1.05; end
    if ~exist('rs','var'); rs = 4.2e-10; end
    if ~exist('T', 'var'); T = 298.15; end
    if ~exist('tau','var'); tau = 5e-2; end
    if ~exist('nsteps','var'); nsteps = 50000; end
    if ~exist('epsilon','var'); epsilon = 1e-6; end
    if ~exist('bsair','var'); bsair = false; end
    if ~exist('scale_LP','var'); scale_LP = true; end
    if ~exist('k0','var'); k0 = 50; end
    if ~exist('writemat','var'); writemat = false; end
    if ~exist('note','var'); note = ''; end
    if ~exist('rdt','var'); rdt = 10; end

    % Set transpiration?
    Eset = exist('E','var');
    if ~exist('EperBS','var'); EperBS = true; end

    % Mesophyll concentrations
    if ~exist('cm', 'var')
        if exist('cm_star','var')
            cm = cm_star;
        else
            cm = 500;
            cm_star = cm;
        end
    elseif ~exist('cm_star','var')
        if length(cm)==1 % If just one value given, all ms set to this value
            if cm==0
                cm_star = 500;
            else
                cm_star = cm;
            end
        else
            cm_star = mean(cm);
        end
    end
    
    % If sigma is not given, calculate using the hindrance factor
    if ~exist('sig','var')
        W = @(lam) (1 - 3.02*lam.^2 + 5.776*lam.^3 - 12.3674*lam.^4 ...
            + 18.9775*lam.^5 - 15.2185*lam.^6 + 4.8525*lam.^7);
        sig = 1 - W(rs/h);
    end
    if ~exist('psig','var'); psig = false; end
    
    % Check if the effective diffusion (omega) is given
    set_omega = exist('omega', 'var');
    
    % Characteristic pressure
    p0 = eta*D/h^2;
    if ~exist('Pb','var'); Pb = [0;1;-15]; end
    Pb = Pb*1e5/p0;
    Pb = reshape(Pb, 3, 1);
    
    % Number of nodes
    N = length(IN);
    IN = reshape(IN,N,1);
    
    % Concentrations
    % Boundary concentrations
    if ~exist('cb','var'); cb = [0; 0]; else
        switch length(cb)
            case 1
                % If only one boundary concentration given, assign to phloem
                cb = [0 cb]';
            case 2
                % Else, full assign
                cb = cb(:);
            otherwise
                error('Invalid cb construction. Must be of length 1 or 2')
        end
        
        % If the phloem concentration exceeds 1 (ie, is "higher than the
        % mesophyll", interpret this as a real concentration instead of a scaled
        % concentration. Rescale by dividing by the mesophyll concentration.
        if cb(2) > 1
            cb = cb/cm_star;
            warning(paste('Boundary phloem concentration interpreted as a',...
                'real concentration. Rescaling by the mesophyll concentration'))
        end
    end
    
    % Permeabilities etc
    % Hindrance factor (movement of sugar through PD)
    H = @(lam) (1 + (9/16).*lam.*log(lam) - 1.19358.*lam + 0.4285.*lam.^3 ...
        - 0.3192.*lam.^4 + 0.08428.^lam.^5);
    H0 = H(rs/h);
    
    % Unpack permeability input structure/vector
    if exist('LP','var')
        % If given as a structure with named fields, fill in vector
        % If individual connections are named, fill in separately and
        % replace missing values with the defaults
        if isstruct(LP)            
            % General Lp vector within the nodal network (see table in notes)
            fields = {'axax','axtt','tttt',...
                'axap','axtp','ttap','tttp','axbs','ttbs',...
                'apap','aptp','tptp','apbs','tpbs','bsbs', 'bsms'};
            
            % Default LPs corresponding to fields above
            default_LP = [1 1 1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1];
            
            LPtemp1 = zeros(16,1);
            LPtemp2 = zeros(16,1);

            % Tracheary connections
            if isfield(LP,'trach'); LPtemp1(1:3) = LP.trach;
            else; LPtemp2(1:3) = default_LP(1:3); end
            
            % Osmotic connections
            if isfield(LP,'osm'); LPtemp1(4:9) = LP.osm;
            else; LPtemp2(4:9) = default_LP(4:9); end

            % Phloem connections
            if isfield(LP,'phloem'); LPtemp1(10) = LP.phloem;
            else; LPtemp2(10) = default_LP(10); end

            % Plasmodesmal connections
            if isfield(LP,'pd'); LPtemp1(11:16) = LP.pd;
            else; LPtemp2(11:16) = default_LP(11:16); end

            % Individually determined permeabilities/conductivities
            for i = 1:16
                if isfield(LP,fields{i})
                    LPtemp1(i) = LP.(fields{i});
                    LPtemp2(i) = 0; 
                end
            end
            
            % Combine default and given values
            LP = LPtemp1 + LPtemp2*max([LPtemp1(:); LPtemp2(:)]);
        else
            % Processing of vector inputs
            if length(LP) == 1
                LP = ones(16,1)*LP;
            elseif length(LP) ~= 16
                clear LP; 
                warning('Invalid format for LP input. Using default values.')
            else
                LP = reshape(LP, 16,1); 
            end
        end
    end
    
    % Default permeabilities/conductivities
    if ~exist('LP','var')
        LP = zeros(16,1);
        LP(1:3) = 1; % Tracheary, phloem
        LP(4:9) = 0.1; % Osmotic
        LP(10:16) = 0.1; % Living/PD
    end
    
    % Then move everyone up by one so that LP can be applied directly to Chi
    LP = [0;LP];
    
    if ~exist('Lp_bounds','var')
        Lp_bounds = [LP(2)*.2;LP(11)*.2;max(LP)*1.5e-3]; 
    else
        Lp_bounds = reshape(Lp_bounds,3,1); 
    end
    
    % If the permeabilities should be scaled
    if scale_LP
        lpscale = max(max(LP),max(Lp_bounds));
        Lp_bounds = Lp_bounds/lpscale;
        LP = LP/lpscale;
        if exist('omega', 'var'); omega = omega/lpscale; end
    else
        lpscale = 1;
    end
    
    % Starch kinetics
    ks_def = struct('vmax',1,'km',1,'kh',1);
    if ~exist('ks', 'var')
        ks = ks_def;
    elseif ~isstruct(ks) && length(ks)==1
        ks = struct('vmax',ks,'km',ks,'kh',ks);
        % If no starch kinetics (ie, all zero), set km to a nonzero value to
        % prevent divisions by zero
        if ks.vmax == 0
            ks.km = 1;
        end
    elseif isstruct(ks)
        ks = renamefield(ks,lower(fieldnames(ks)));
        for i = {'vmax','kh','km'}
            if ~ismember(i{1},fieldnames(ks))
                ks.(i{1}) = ks_def.(i{1});
            end
        end
    else
        warning('Invalid starch parameters. Using default values.')
        ks = ks_def;
    end

    %% General setup
    
    % Extract original connection matrix
    chi = triu(chi);
    
    % Find "edge" xylem and phloem (ie, axial xylem and phloem that are actually
    % connected to the boundary nodes
    % Xylem
    is_ex = strcmp(IN,'ex');
    if sum(is_ex)~=0
        IN{is_ex} = 'ax';
    else
        is_ex = strcmp(IN, 'ax');
    end
    
    % Phloem
    is_ep = strcmp(IN,'ep');
    if sum(is_ep)~=0
        IN{is_ep} = 'ap';
    else
        is_ep = strcmp(IN, 'ap');
    end
    
    % Convert id strings to prime indices
    primes = [2 3 5 7 11 13];
    ID = zeros(length(IN),1);
    ID(strcmp('ax',IN)) = primes(1);
    ID(strcmp('tt',IN)) = primes(2);
    ID(strcmp('ap',IN)) = primes(3);
    ID(strcmp('tp',IN)) = primes(4);
    ID(strcmp('bs',IN)) = primes(5);
    
    % Reorder in case the input vector wasn't properly ordered
    [ID,ind] = sort(ID);    
    chi = chi(ind,ind);
    IN = IN(ind);
    
    % Initial concentrations
    if ~exist('c0','var')
        c0 = zeros(N,1);
    else
        if (length(c0)==N)
            c0 = reshape(c0, N, 1);
        else
            nt = sum(ismember(ID,[2, 3]));
            if (length(c0)==1) || (length(c0)==(N-nt))
                temp = zeros(N,1);
                temp((nt+1):end) = c0;
                c0 = temp;
            else
                warning('Invalid c0 construction. Setting all initial concentrations to zero.')
            end
        end
    end
    
    % Initial starch
    if ~exist('s0','var')
        s0 = zeros(N,1);
    else
        if (length(s0)==N)
            s0 = reshape(s0,N,1);
        else
            n_nostarch = sum(ismember(ID,[2,3,5]));
            if (length(s0)==1) || (length(s0)==(N-n_nostarch))
                temp = zeros(N,1);
                temp((n_nostarch+1):end) = s0;
                s0 = temp;
            else
                warning('Invalid s0 construction. Setting all initial starch concentrations to zero.')
            end
        end
    end
    
    % Find number of bundle sheath cells
    nbs = sum(ID==primes(5));
    
    % Add in mesophyll cells
    ID((end+1):(end+nbs)) = primes(6);
    is_ex((end+1):(end+nbs)) = 0;
    is_ep((end+1):(end+nbs)) = 0;

    % Scale ET if given as E per BS cell
    if Eset && EperBS
        E = E*nbs;
    end

    
    % Mesophyll concentration
    if length(cm) == 1
        if cm == 0
            cm = zeros(nbs,1);
        else
            cm = ones(nbs,1);
        end
    elseif length(cm) == nbs
        if abs(mean(cm)-cm_star)/cm_star <= 0.1
            cm = reshape(cm,nbs,1)/cm_star;
        else
            cm = reshape(cm,nbs,1);
        end
    else
        warning('Mismatched cm. Setting all mesophyll concentrations to 1.')
        cm = ones(nbs,1);
    end
    cm0 = cm;
    
    % Expand chi to include bs-m connections
    chi((N-nbs+1):N,(N+1):(N+nbs)) = eye(nbs);
    % Fill in bottom rows to make chi square
    chi((N+1):(N+nbs),:) = 0;
    % Reflect chi to get symmetry
    chi = chi + chi';

    % Get unique codes for each connection type
    Chi = chi.*ID.*ID';

    % Collapse into simple vector (cvt for "convert")
    cvt = sparse(169,1);
    cvt([unique(primes(1:5).*primes(1:5)'); primes(5)*primes(6)]) = [1:8 10 9 11:16];
    Chi(Chi>0) = cvt(Chi(Chi>0)); 

    % Get the Lp for each connection based on connection type
    xi = LP(Chi+1).*chi;

    % Boundary conditions
    % Determine boundary-adjacent nodes (ax, ap, bs)
    atbounds = struct();
    atbounds.bool = ismember(ID,[11 13]) | is_ex | is_ep;
    
    % atbounds.which: boolean matrix that indicates which boundary is adjacent
    % to a given node-- whichbounds * cb is an Nx1 vector giving the
    % concentration of the boundary corresponding to each node (0 if
    % internal)
    atbounds.which = zeros(N+nbs,2);
    atbounds.which(is_ex,1) = 1; % Xylem
    atbounds.which(is_ep,2) = 1; % Phloem
    
    % atbounds.w is a (N+nbs)x1 vector that gives the permeability between a 
    % node and its corresponding boundary. Bundle sheath (ID==11) is only given 
    % the permeability to the air space, since the connection to the mesophyll 
    % is now explicitly modeled. Similarly, the mesophyll (ID==13) is only given 
    % the permeability to the air space.
    atbounds.w = zeros(N+nbs,1);
    atbounds.w(is_ex) = Lp_bounds(1);
    atbounds.w(is_ep) = Lp_bounds(2);
    atbounds.w(ID==11) = Lp_bounds(3) * (1-bsair);
    atbounds.w(ID==13) = Lp_bounds(3);
    
    % atbounds.omega (Nx1) gives the effective diffusion across an interface at
    % the boundaries. For the xylem and phloem, this is just 3*H0*w. For
    % the bundle sheath, since there is no sugar flow into the airspace, this is
    % just zero. Mesophyll is ignored here.
    if set_omega
        atbounds.omega(ismember(ID(1:N), [11 13]) | is_ex(1:N) | is_ep(1:N),1) = omega;
    else
        atbounds.omega = 3*H0*atbounds.w(1:N);
    end
    atbounds.omega(ID==11) = 0;

    % Make Umat (internal pressure connections)
    % Set identity diagonal to balance out other nodes, then also subtract
    % contributions from boundary elements
    Umat = (xi - (sum(xi,2)+atbounds.w).*eye(N+nbs));
    if Eset
        Umat(:,end+1) = atbounds.w.*(ismember(ID,[11 13]));
        Umat(end+1,1:(end-1)) = Umat(:,end)';
        Umat(end,end) = -sum(Umat(end,:));
    end


    % Define sigmas (reflection coefficients)
    % sigma = 1 means impermeability, so big osmotic effect
    % sigma = 0 means complete permeability, so no osmotic effect
    sigma = struct('trach',0, 'osm',1, 'pd', sig, ...
        'xylem', 0, 'phloem', sig*psig, 'meso', sig);
    
    % Main Sigma matrix assigns sigma depending on connection type WITHIN
    % network (ie, no boundaries), so mesophyll boundary is PD w/ bundle sheath
    Sigma = Chi*0;
    Sigma(ismember(Chi, 1:3)) = sigma.trach;    % tracheary
    Sigma(ismember(Chi, 4:9)) = sigma.osm;      % osmotic
    Sigma(ismember(Chi, 10:15)) = sigma.pd;     % plasmodesmal
    Sigma(Chi==16) = sigma.meso;                % mesophyll
    
    %% Boundary conditions/constraints
    
    % Boundary pressure constraints:
    % For xylem and phloem, this should be permeability * pressure at the
    % boundary.
    % For bundle sheath, this will only include the air space since the
    % mesophyll pressure is modeled explicitly.
    % Same for mesophyll; this constraint is only with the air space.
    u_p = zeros(N+nbs,1);
    u_p(is_ex) = Lp_bounds(1)*Pb(1);            % Xylem
    u_p(is_ep) = Lp_bounds(2)*Pb(2);            % Phloem
    if ~Eset % If E is not set, use constant air pressure
        u_p(ID==11) = Lp_bounds(3)*Pb(3)*(1-bsair); % Bundle sheath
        u_p(ID==13) = Lp_bounds(3)*Pb(3);           % Mesophyll
    else
        u_p(end+1) = -abs(E);
    end
    
    % Sigma values at boundary-adjacent nodes
    % (not) Multiplied by their respective permeabilities
    atbounds.Sigma = atbounds.which;
    atbounds.Sigma(is_ex,1) = sigma.xylem;%*Lp_bounds(1);
    atbounds.Sigma(is_ep,2) = sigma.phloem;%*Lp_bounds(2);
    atbounds.Sigma(ID==11,3) = sigma.osm;%*Lp_bounds(3)*(1-bsair);  % with air
    atbounds.Sigma(ID==13,3) = sigma.osm;%*Lp_bounds(3);     % with air
    
    % A condensed form of sigma for use with sugar flow calculations (ie,
    % each node is adjacent to at most one boundary)
    atbounds.sigma = sum(atbounds.Sigma(1:N,:),2);
    
    atbounds.Sigma = atbounds.Sigma .* atbounds.w;
    
    % Find plasmodesmal connections
    non_osm = ismember(Chi, [1:3 10:16]);
    
    % Define Cmat: for INTERNAL connections only
    Cmat = xi.*Sigma - eye(N+nbs).*(sum(xi.*Sigma,2) + sum(atbounds.Sigma,2));

    % Add last row if transpiration BC is set
    if Eset
        Cmat(end+1,:) = sum(atbounds.Sigma,2)'.*(ismember(ID,[11 13]))';
        Cmat(:,end+1) = 0; % Since c_air = 0
        atbounds.Sigma(end+1,:) = 0;
    end
    
    % Remove last column of atbounds.Sigma since it's no longer necessary
    atbounds.Sigma = atbounds.Sigma(:,1:2);
    
    % Extract the mesophyll columns of Cmat because cm is not included in c
    cm_corr = Cmat(:,(N+1):(N+nbs));
    
    % Remove last columns of Cmat (redundant because of constant cm)
    Cmat = Cmat(:,1:N);

    % Define wRT for sugar flux calculations
    RT = 8.314*T*cm_star*h^2/(eta*D);
    if set_omega
        wRT = omega.*non_osm.*RT;
    else
        wRT = 3*H0*xi.*non_osm.*RT;
    end
    
    % Mesophyll boundary condition
    is_ms = (N+1):(N+nbs);
    if ~exist('jm','var') || isnan(jm)
        if allow_buildup
            cm_update = @(c_m, J_ij, dt)...
                (min(maxcm*cm0,max(cm0, c_m + sum(J_ij(is_ms,:),2)*dt)));
        else
            cm_update = @(c_m, J_ij, dt) (c_m);
        end
    else
        cm_update = @(c_m, J_ij, dt) (c_m + (jm + sum(J_ij(is_ms,:),2))*dt);
    end
    
    %% Initializing structures and values
    
    % Initialize P and C
    % For in-network nodes
    P = zeros(ceil(nsteps/rdt),N);
    C = zeros(ceil(nsteps/rdt),N);
    S = zeros(ceil(nsteps/rdt),N);
    
    % Concentration of the boundary phloem and mesophyll
    Cp = zeros(ceil(nsteps/rdt), 1);
    Cm = zeros(ceil(nsteps/rdt), nbs);

    % Sugar export
    Jp = zeros(ceil(nsteps/rdt),1);
    
    % Time
    time = 0;
    Time = zeros(ceil(nsteps/rdt),1);
    
    % Pressure of the mesophyll (and air, if E set)
    Pm = zeros(ceil(nsteps/rdt), nbs);
    if Eset
        Pa = zeros(ceil(nsteps/rdt),1);
    end

    % Initialize loop
    c = c0; s = s0;
    l = 0; % step counter
    k = 0; % readout counter
    loop = true; % boolean that updates w/ convergence or max steps
    opts.SYM = true; %for linsolve
    
    % cell volume %
    vcell = 1;
    
    % Which cells to calculate starch kinetics for
    for_starch = strcmp(IN,'tp') | strcmp(IN,'bs');
    
    % Shorthand for the main step
    mainstep = @(c,cm,s,cb) systemsolve(c,cm,s,...
            chi,N,RT,Cmat,cm_corr,atbounds,cb,u_p,Umat,opts,Sigma,xi,wRT,for_starch,vcell,ks);
        
    % Time epsilon
    teps = 2*epsilon;
    if epsilon==0
        teps = 2e-6;
    end
    
    %% Main loop
    while loop
        % Increment k
        l = l+1;
        
        %-----------------------------------------------------------------------
        % Adaptive time step
        % Calculate values at tk
        [p,Uij,Ucij,wRTc,Jij,Uib,Jib,Ji,dsdt] = mainstep(c,cm,s,cb);
        
        % One full time step
        c1 = c + (Ji-dsdt)*tau;
        cm1 = cm_update(cm,Jij,tau);
        s1 = s + dsdt*tau;
        
        % Half time step
        c05 = c + (Ji-dsdt)*tau/2;
        cm05 = cm_update(cm,Jij,tau/2);
        s05 = s + dsdt*tau/2;
        cp05 = 0;
        
        % Second half time step
        [~,~,~,~,Jij2,~,~,Ji2,dsdt2] = mainstep(c05,cm05,s05,[cb(1);cp05]);
        c2 = c05 + (Ji2-dsdt2)*tau/2;
        cm2 = cm_update(cm05,Jij2,tau/2);
        s2 = s05 + dsdt2*tau/2;
        
        % Determine time step factor
        err = [c2;cm2;s2] - [c1;cm1;s1];%[c2;cm2;s2;cp2] - [c1;cm1;s1;cp1];
        if norm(err(1:(end-1))) <= 5*teps
            tfac = .95*min(max((teps/(norm(err(end))))^0.5,0.3),2);
        else
            tfac = .95*min(max((teps/(norm(err)))^0.5,0.3),2);
        end
        
        % Update the time step
        tau = tau*tfac;
        
        % Do the sugar update with the new time step
        c = c + (Ji-dsdt)*tau;
        cm = cm_update(cm,Jij,tau);
        s = s + dsdt*tau;
        
        % Add time step to the total time counter
        time = time + tau;
        
        %-----------------------------------------------------------------------
        % Save information
        if mod(l,rdt)==0
            k = k+1;
            Cp(k) = mean(c(is_ep));
            P(k,:) = p(1:N)';
            Pm(k,:) = p((N+1):(N+nbs))';
            if Eset; Pa(k) = p(end); end
            C(k,:) = c';
            S(k,:) = s';
            Cm(k,:) = cm';
            Jp(k) = sum(Jib(is_ep));
            Time(k) = time;
        end
        
        %-----------------------------------------------------------------------        
        % Check for convergence
        if epsilon == 0
        elseif any(c<0) || (k>k0 && (norm(Ji-dsdt) <= epsilon) && (abs((Cp(k)-Cp(k-k0+1))/tau) <= epsilon))
            loop = false;
            
            % Truncate if converged early
            C = C(1:k,:);
            P = P(1:k,:);
            S = S(1:k,:);
            Cp = Cp(1:k);
            Pm = Pm(1:k,:);
            if Eset; Pa = Pa(1:k); end
            Cm = Cm(1:k,:);
            Jp = Jp(1:k);
            Time = Time(1:k);

            if any(c<0)
                warning('Negative concentrations. Ending simulation.');
                C(end) = NaN;
                P(end) = NaN;
                S(end) = NaN;
                Cp(end) = NaN;
                Pm(end) = NaN;
                Cm(end) = NaN;
                Uij = Uij*NaN;
                Jij = Jij*NaN;
                Uib = Uib*NaN;
                Jib = Uib*NaN;
            end
        end
        
        if (l>=nsteps)
            loop = false;
        end
    end % Main loop end
    
    % Realign Cp
%     Cp = [Cp(2:k); cb(2)];
    
    %% Write node and edge tables to file if filename is given
    if exist('filename','var')
        fprintf('Simulation finished. Writing data...\n')
        
        %-----------------------------------------------------------------------
        % Nodes
        % Pre-assigning node values
        ms = mat2cell(repelem('ms',nbs,1),nbs,2);
        names = cat(1,{'bx';'bp'},IN,cellstr(ms{1}));
        pressures = [Pb(1:2); p(1:(N+nbs))];
        concentrations = [cb; C(end,:)'; cm];
        starches = [0; 0; S(end,:)'; zeros(nbs,1)];
        
        % Build node table
        nodeOp = table((-1:(N+nbs))', names, ...
            pressures, concentrations, starches, ...
            RT*concentrations,...
            pressures - RT*concentrations,...
            'VariableNames', {'Id','Label', 'P', 'C', 'S', 'Pi', 'Psi'});
        
        % Additional nodes for unmodeled airspace
        esnodes = repmat(cell2table({-200, 'a', Pb(3), 0, 0, 0, Pb(3)},...
            'VariableNames',{'Id','Label', 'P', 'C', 'S', 'Pi', 'Psi'}),nbs,1);
        esnodes.Id = esnodes.Id-find(ID==11);
        
        % Assemble full node table
        nodeOp = [nodeOp;esnodes];
        
        %-----------------------------------------------------------------------
        % Edges
        % Edge table column names
        varnames = {'Source','Target','flow','value'};
        
        % Pre-allocate space for edge table
        nedge_start = 2*sum(chi(:)) + sum(ismember(ID,[2 5 11 13]));
        edgeOp = cell2table(repmat({0, 0, 'x', 0},nedge_start,1), 'VariableNames', varnames);
        l = 1;
        
        for i = 1:(N+nbs)
            % Go through each internal connection
            for j = (i+1):(N+nbs)
                if (Chi(i,j)==0)
                    continue
                    % ignore if not actually connected
                end
                
                % Make sure the correct direction is specified
                % Water flux
                if Uij(i,j) > 1e-6
                    edgeOp(l,:) = cell2table({j,i,'u',Uij(i,j)}, 'VariableNames', varnames);
                    l = l + 1;
                elseif Uij(i,j) < -1e-6
                    edgeOp(l,:) = cell2table({i,j,'u',-Uij(i,j)}, 'VariableNames', varnames);
                    l = l + 1;
                end
                % Sugar (or solute in general) flux
                if Jij(i,j) > 1e-6
                    edgeOp(l,:)  = cell2table({j,i,'j',Jij(i,j)}, 'VariableNames', varnames);
                    l = l + 1;
                elseif Jij(i,j) < -1e-6
                    edgeOp(l,:)  = cell2table({i,j,'j',-Jij(i,j)}, 'VariableNames', varnames);
                    l = l + 1;
                end
            end
            
            % Edges from boundary nodes
            switch ID(i)
                case {2,5} % Xylem and phloem, but only edge nodes ex & ep
                    if (is_ex(i) || is_ep(i))
                        bnode = (ID(i)-5)/3;
                        if Uib(i) > 1e-6
                            edgeOp(l,:)  = cell2table({bnode,i,'u',Uib(i)},'VariableNames',varnames);
                            l = l + 1;
                        elseif Uib(i) < -1e-6
                            edgeOp(l,:)  = cell2table({i,bnode,'u',-Uib(i)},'VariableNames',varnames);
                            l = l + 1;
                        end
                        if Jib(i) > 1e-6
                            edgeOp(l,:)  = cell2table({bnode,i,'j',Jib(i)},'VariableNames',varnames);
                            l = l + 1;
                        elseif Jib(i) < -1e-6
                            edgeOp(l,:)  = cell2table({i,bnode,'j',-Jib(i)},'VariableNames',varnames);
                            l = l + 1;
                        end
                    end
                case {11,13} % Bundle sheath and mesophyll
                    % Each bundle sheath shares one (visual) airspace node with
                    % its respective mesophyll cell
                    if ID(i) == 11
                        bnode = -200-i;
                    else
                        bnode = -200-(i-nbs);
                    end
                    if Uib(i) > 1e-6
                        edgeOp(l,:)  = cell2table({bnode,i,'u',Uib(i)},'VariableNames',varnames);
                        l = l + 1;
                    elseif Uib(i) < -1e-6
                        edgeOp(l,:)  = cell2table({i,bnode,'u',-Uib(i)},'VariableNames',varnames);
                        l = l + 1;
                    end
            end
        end
        
        % Remove extra rows
        edgeOp = edgeOp(1:(l-1),:);
        
    
        % Reassign edge identities to ex and ep
        IN{is_ex} = 'ex';
        IN{is_ep} = 'ep';
        
        %-----------------------------------------------------------------------
        % Write tables to files
        writetable(nodeOp, strjoin({'output/',filename,'_nodes.csv'},''))
        writetable(edgeOp, strjoin({'output/',filename,'_edges.csv'},''))
        fprintf('Data written to \n%s\n%s\n', ...
            strjoin({filename,'_nodes.csv'},''), ...
            strjoin({filename,'_edges.csv'},''))
    end
    
    % Reassign edge identities to ex and ep
    IN{is_ex} = 'ex';
    IN{is_ep} = 'ep';
    
    %% Create output structure
    
    output = struct('C',C,'P',P,'S',S,'Time',Time,'Psi', P - C*RT, ...
        'Jp',Jp,'Uij',Uij,'Jij',Jij,'Uc',Ucij,'Dc',wRTc,...
        'Uib',Uib,'Jib',Jib,...
        'k',k,'Cp',Cp, 'Pm', Pm, 'Cm', Cm, 'ks', ks, ...
        'RT',RT, 'sigma', sig, 'IN', {IN}, 'chi', chi, 'Chi', Chi, ...
        'Lp', LP(2:end), 'Lp_bounds', Lp_bounds, 'lpscale', lpscale, ...
        'note',note,'args', args,...
        'SIconvert', @(toconvert)(toSI(toconvert,h,eta,D,cm_star)));
    if Eset; output.Pa = Pa; end
    
    %% Save output structure if filename given
    if exist('filename', 'var') && writemat
        runname = strjoin({filename,'_run'},'');
        eval([runname '= output;']);
        save(strjoin({'output/',filename,'.mat'},''), runname)
        fprintf('%s\n',strjoin({filename,'.mat'},''))
    end
end

function [p,Uij,Ucij,wRTc,Jij,Uib,Jib,Ji,dsdt] = systemsolve(c,cm,s,...
    chi,N,RT,Cmat,cm_corr,atbounds,cb,u_p,Umat,opts,Sigma,xi,wRT,for_starch,vcell,ks)

    Ck = [c;cm];
    Ckt = [c;cm];
    %-----------------------------------------------------------------------
    % Pressure calculation
    % Get osmotic constraints
    u_osm = RT*(Cmat*c + cm_corr*cm + atbounds.Sigma*cb);

    % Add boundary pressure and osmotic constraints
    b = -(u_p - u_osm);

    % Solve for p
    p = linsolve(Umat, b, opts);

    if length(p) > length(Ckt)
        pa = p(end);
        p = p(1:(end-1));
    end

    %-----------------------------------------------------------------------
    % Calculate individual fluxes
    %---Advective flux---
    % Advective speed
    Uij = xi.*(p'-p - Sigma.*RT.*(Ckt'-Ckt));

    if exist('pa','var')
        p = [p; pa];
    end

    % Multiply by pairwise mean concentration to get sugar flux
    Ucij = 0.5*Uij.*(Ck'+Ck).*(1-Sigma).*chi;

    %---Diffusive fluxes---
    wRTc = wRT.*(Ck'-Ck);

    % Add together
    Jij = Ucij + wRTc;

    %---Flow to boundaries---
    % Use water balance to find the flow going to/from boundaries
    Uib = -sum(Uij.*(atbounds.bool),2);

    % Determine the relevant boundary concentration
    c_use = 0.5*(atbounds.which(1:N,:)*cb + c);

    % Calculate sugar flux at the boundaries
    Jib = ((1-atbounds.sigma).*Uib(1:N)).*c_use;% + ...
       % atbounds.omega.*(atbounds.which(1:N,:)*cb - c)*RT);
    
    %-----------------------------------------------------------------------
    % Sugar update
    % Calculate sugar fluxes
    Ji = sum(Jij(1:N,:),2) + Jib;
    
    % Starch kinetics
    dsdt = for_starch.*(vcell.*(ks.vmax.*c)./(ks.km + c) - ks.kh.*s);
end
