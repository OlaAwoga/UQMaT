classdef Params
    % UQMaT.Params  All system parameters in one place.
    %
    %   p = UQMaT.Params('system', <name>, Name, Value, ...)
    %
    % ====================================================================
    %  SYSTEM
    % ====================================================================
    %   'chain'   : 1D uniform tight-binding chain
    %   'ssh'     : 1D Su-Schrieffer-Heeger chain
    %   'kitaev'  : 1D Kitaev p-wave superconductor
    %   'square'  : 2D square lattice
    %
    % ====================================================================
    %  HAMILTONIAN SIZE  (set automatically)
    % ====================================================================
    %   sc=false, spin=false  →  N × N
    %   sc=false, spin=true   →  2N × 2N
    %   sc=true               →  4N × 4N   BdG: Ψ=(c↑,c↓,c↑†,c↓†)ᵀ
    %   kitaev, Delta_p≠0     →  2N × 2N   spinless BdG
    %   kitaev, Delta_p=0     →  N × N
    %
    % ====================================================================
    %  GEOMETRY
    % ====================================================================
    %   N        : sites (unit cells for ssh)              (default: 50)
    %   bc       : 'obc' | 'pbc'                           (default: 'obc')
    %   Nx, Ny   : sites along x, y  (square)              (default: 10)
    %   bcx, bcy : 'obc'|'pbc' per direction  (square)     (default: 'obc')
    %
    % ====================================================================
    %  HOPPING
    % ====================================================================
    %   t        : hopping amplitude                        (default: 1.0)
    %              scalar    → uniform (all lattices)
    %              [tx, ty]  → anisotropic (square lattice only)
    %   mu       : uniform chemical potential               (default: 0.0)
    %   onsite   : uniform on-site energy; H(i,i) = onsite - mu (default: 0.0)
    %   t2       : next-nearest hopping                     (default: 0.0)
    %   v        : SSH intra-cell hopping                   (default: 0.5)
    %   w        : SSH inter-cell hopping                   (default: 1.0)
    %
    % ====================================================================
    %  ZEEMAN + RASHBA
    % ====================================================================
    %   B        : Zeeman field [Bx, By, Bz]               (default: [0,0,0])
    %   lambda_R : Rashba SOC strength                      (default: 0.0)
    %
    % ====================================================================
    %  SUPERCONDUCTIVITY
    % ====================================================================
    %   sc       : enable BdG (auto-set when any Delta nonzero) (default: false)
    %   Delta_s  : singlet s-wave amplitude                 (default: 0.0)
    %              scalar       → uniform
    %              N-vector     → site-resolved (1D)
    %              Nx×Ny matrix → spatial map (2D)
    %   dvec     : triplet d-vector [dx, dy, dz]            (default: [0,0,0])
    %   Delta_t  : triplet bond amplitude                    (default: 0.0)
    %   Delta_p  : Kitaev p-wave amplitude                   (default: 0.0)
    %              scalar | (N-1)-vector
    %   pairing_mode is auto-inferred: 'singlet'|'triplet'|'mixed'
    %
    % ====================================================================
    %  DISORDER
    % ====================================================================
    %   W        : disorder strength (W>0 auto-enables 'anderson') (default: 0.0)
    %   disorder : disorder type                             (default: 'none')
    %     'anderson'    : ε_i ~ Uniform[-W/2,+W/2], zero mean
    %     'gaussian'    : ε_i ~ N(0,W²), zero mean
    %     'binary'      : ε_i ∈ {-W/2,+W/2} (bimodal alloy)
    %     'hop_disorder': δt_ij random on all bonds (off-diagonal)
    %     'impurity'    : controlled impurities (see imp_* params)
    %     'vacancy'     : sites at imp_pos removed
    %   seed     : RNG seed (0 = random each run)           (default: 0)
    %
    %   Impurity parameters:
    %   imp_pos      : position(s): 1D site index/vec; 2D [ix,iy] or Nimp×2
    %   imp_potential: on-site potential V₀ (scalar or Nimp-vector)
    %   imp_magnetic : exchange strength J  (scalar or Nimp-vector)
    %   imp_sigma    : Gaussian spread σ in lattice units (0 = point)
    %   imp_geometry : 2D shape: 'point'|'hline'|'vline'|'line'|'cluster'
    %
    % ====================================================================
    %  SPECTRAL
    % ====================================================================
    %   eta        : Lorentzian broadening                  (default: 0.01)
    %   EnergyGrid : energy grid vector, e.g. -2:0.005:2   (default: auto)
    %
    % ====================================================================
    %  EXAMPLES
    % ====================================================================
    %   % Normal chain
    %   p = UQMaT.Params('system','chain','N',80,'t',1,'mu',0)
    %
    %   % Rashba + Zeeman chain (spinful, 2N×2N)
    %   p = UQMaT.Params('system','chain','N',80,'B',[0,0,0.3],'lambda_R',0.5)
    %
    %   % s-wave SC chain (BdG, 4N×4N, sc auto-enabled)
    %   p = UQMaT.Params('system','chain','N',80,'Delta_s',0.3)
    %
    %   % Topological SC candidate (Rashba + Zeeman + Delta_s)
    %   p = UQMaT.Params('system','chain','N',80, ...
    %       'Delta_s',0.3,'B',[0,0,0.5],'lambda_R',0.6)
    %
    %   % SSH
    %   p = UQMaT.Params('system','ssh','N',40,'v',0.3,'w',1.0)
    %
    %   % Kitaev (BdG when Delta_p ≠ 0)
    %   p = UQMaT.Params('system','kitaev','N',60,'mu',-1,'Delta_p',0.5)
    %
    %   % 2D square lattice
    %   p = UQMaT.Params('system','square','Nx',20,'Ny',20)
    %
    %   % 2D square + SC + anisotropic hopping
    %   p = UQMaT.Params('system','square','Nx',16,'Ny',16, ...
    %       'Delta_s',0.3,'t',[1.0,0.5])
    %
    %   % Custom energy grid
    %   p = UQMaT.Params('system','chain','N',60,'Delta_s',0.3, ...
    %       'EnergyGrid',-2:0.005:2)
    %
    %   % Site-resolved Delta_s (proximity effect)
    %   Ds = zeros(60,1); Ds(20:40) = 0.3;
    %   p = UQMaT.Params('system','chain','N',60,'Delta_s',Ds)
    %
    %   % Anderson disorder
    %   p = UQMaT.Params('system','chain','N',100,'W',2.0,'seed',42)
    %
    %   % Point impurity at site 30
    %   p = UQMaT.Params('system','chain','N',80,'Delta_s',0.3, ...
    %       'imp_pos',30,'imp_magnetic',2.0)
    %
    %   % 2D horizontal line impurity
    %   p = UQMaT.Params('system','square','Nx',20,'Ny',20, ...
    %       'imp_pos',[1,10],'imp_potential',1.0,'imp_geometry','hline')

    %  Citation
    %  --------
    %  UQMaT: see CITATION.cff in the package root.
    %  Tight-binding method: Datta S., "Electronic Transport in Mesoscopic
    %    Systems", Cambridge (1995).
    %  BdG formalism: Zhu J.-X., "Bogoliubov-de Gennes Method and Its
    %    Applications", Springer (2016).
    properties
        system   (1,:) char  = 'chain'

        % Geometry 1D
        N        (1,1) double {mustBeInteger,mustBePositive} = 50
        bc       (1,:) char  = 'obc'

        % Geometry 2D
        Nx       (1,1) double {mustBeInteger,mustBePositive} = 10
        Ny       (1,1) double {mustBeInteger,mustBePositive} = 10
        bcx      (1,:) char  = 'obc'
        bcy      (1,:) char  = 'obc'

        % Geometry 3D
        Nz       (1,1) double {mustBeInteger,mustBePositive} = 10
        bcz      (1,:) char  = 'obc'

        % Heterostructure
        hetero   struct = struct()   % heterostructure definition (see Heterostructure classes)

        % TMD material (for system='tmd')
        tmd_material (1,:) char = 'MoS2'   % 'MoS2'|'WS2'|'MoSe2'|'WSe2'|'custom'
        tmd_params   struct = struct()      % custom TMD parameters

        % Quasiclassical (Eilenberger/Usadel/Keldysh)
        qc_geometry (1,:) char   = 'bulk'  % 'bulk'|'slab'|'ns_junction'|'sns_junction'|'vortex'
        qc_D        (1,1) double = 1.0     % diffusion constant (Usadel)
        qc_vF       (1,1) double = 1.0     % Fermi velocity (Eilenberger)
        qc_Delta0   (1,1) double = 0.1     % initial/bulk gap
        qc_lambda   (1,1) double = 1.0     % pairing interaction
        qc_T        (1,1) double = 0.01    % temperature
        qc_Lx       (1,1) double = 10.0    % system length
        qc_Nx       (1,1) double {mustBeInteger,mustBePositive} = 200
        qc_N_omega  (1,1) double {mustBeInteger,mustBePositive} = 100
        qc_l_b      (1,1) double = Inf     % Kupriyanov-Lukichev barrier length
        qc_V_L      (1,1) double = 0.0    % left voltage (Keldysh)
        qc_V_R      (1,1) double = 0.0    % right voltage (Keldysh)
        qc_sigma_N  (1,1) double = 1.0    % normal conductance (Keldysh)
        qc_formalism (1,:) char  = 'dirty' % 'dirty' (Usadel) | 'clean' (Eilenberger)

        % Hopping
        t        (1,:) double = 1.0   % scalar or [tx,ty]
        mu       (1,1) double = 0.0
        onsite   (1,1) double = 0.0
        t2       (1,1) double = 0.0
        v        (1,1) double = 0.5   % SSH intra-cell
        w        (1,1) double = 1.0   % SSH inter-cell

        % Zeeman / Rashba
        B        (1,3) double = [0 0 0]
        lambda_R (1,1) double = 0.0

        % Non-Hermitian asymmetry (Hatano-Nelson)
        t_nh     (1,1) double = 0.0   % NH hopping: t_R = t - t_nh, t_L = t + t_nh

        % Superconductivity struct
        % SC = struct('Delta_s',0.3,'sym','onsite','Delta_t',0.1,'sym_t','px')
        % Fields: Delta_s, sym (singlet sym), Delta_t, sym_t (triplet sym)
        SC       struct = struct()

        % Altermagnetism struct
        % AM = struct('J',0.5,'sym','gwave')
        % Also supports multiple components:
        % AM = struct('J_s',0.3,'sym_s','onsite','J_dxy',0.2,'sym_dxy','dxy')
        AM       struct = struct()

        % Charge density wave struct
        % CDW = struct('Delta',0.2,'sym','onsite')
        CDW      struct = struct()

        % Kitaev p-wave (stays as scalar — specific to Kitaev chain)
        Delta_p  double        = 0.0

        % Disorder
        disorder     (1,:) char   = 'none'
        W            (1,1) double = 0.0
        seed         (1,1) double {mustBeInteger,mustBeNonnegative} = 0
        imp_pos      double       = []
        imp_potential double      = 0.0
        imp_magnetic  double      = 0.0
        imp_sigma    (1,1) double = 0.0
        imp_geometry (1,:) char   = 'point'

        % Spectral
        eta        (1,1) double {mustBePositive} = 0.01
        EnergyGrid (1,:) double = []
    end

    % Derived — computed in constructor, read-only thereafter
    properties (SetAccess = private)
        spin         (1,1) logical = false
        is_2d        (1,1) logical = false
        pairing_mode (1,:) char    = 'singlet'
        tx           (1,1) double  = 1.0
        ty           (1,1) double  = 1.0
        tz           (1,1) double  = 1.0
        is_3d        (1,1) logical = false
        onsite_val   (1,1) double  = 0.0
        imp_potential_val (:,1) double = 0.0
        imp_magnetic_val  (:,1) double = 0.0
        Delta_s_vec  double
        Delta_p_vec  (:,1) double
        Delta_t_val  (1,1) double = 0.0

        % sc flag — true when any pairing amplitude is nonzero
        % Derived from SC struct or Delta_p; read by all model classes
        sc           (1,1) logical = false

        % Resolved SC fields (from SC struct)
        sc_Delta_s   double       = 0.0   % singlet amplitude (scalar or N-vector)
        sc_sym_s     (1,:) char   = 'onsite'
        sc_Delta_t   (1,1) double = 0.0   % triplet amplitude
        sc_sym_t     (1,:) char   = 'px'

        % Resolved AM components (from AM struct) — cell array of structs
        % Each element: struct('J', amplitude, 'sym', label)
        am_components = {}

        % Resolved CDW fields (from CDW struct)
        cdw_Delta    (1,1) double = 0.0
        cdw_sym      (1,:) char   = 'onsite'
    end

    methods
        function obj = Params(varargin)
            if mod(numel(varargin),2) ~= 0
                error('UQMaT:Params:badInput','Arguments must be Name-Value pairs.');
            end
            for k = 1:2:numel(varargin)
                name = varargin{k};
                val  = varargin{k+1};
                if ~isprop(obj, name)
                    error('UQMaT:Params:unknownParam', ...
                        'Unknown parameter: "%s"\n\nValid parameters:\n%s', ...
                        name, UQMaT.Params.paramList());
                end
                obj.(name) = val;
            end

            % Normalise strings
            obj.system       = lower(strtrim(obj.system));
            obj.bc           = lower(strtrim(obj.bc));
            obj.bcx          = lower(strtrim(obj.bcx));
            obj.bcy          = lower(strtrim(obj.bcy));
            obj.disorder     = lower(strtrim(obj.disorder));
            obj.imp_geometry = lower(strtrim(obj.imp_geometry));
            % Parse SC/AM/CDW structs into resolved private fields
            obj = obj.parseSC();
            obj = obj.parseAM();
            obj = obj.parseCDW();

            % Validate system
            valid_sys = {'chain','ssh','kitaev','hnkitaev','square','kagome',...
                         'cubic','diamond','tmd',...
                         'junction1d','chain1d_square2d'};
            if ~ismember(obj.system, valid_sys)
                error('UQMaT:Params:badSystem', ...
                    'system must be one of: %s', strjoin(valid_sys,', '));
            end

            % 2D flag
            obj.is_2d = ismember(obj.system,{'square','kagome','tmd'});
            obj.is_3d = ismember(obj.system,{'cubic','diamond'});

            % Normalise bcz
            obj.bcz = lower(strtrim(obj.bcz));

            % Resolve tx, ty, tz from t
            if isscalar(obj.t)
                obj.tx = obj.t;  obj.ty = obj.t;  obj.tz = obj.t;
            elseif numel(obj.t) == 2
                obj.tx = obj.t(1);  obj.ty = obj.t(2);  obj.tz = obj.t(1);
            elseif numel(obj.t) == 3
                obj.tx = obj.t(1);  obj.ty = obj.t(2);  obj.tz = obj.t(3);
            else
                error('UQMaT:Params:badT', ...
                    't must be scalar, [tx,ty] or [tx,ty,tz]. Got length %d.', numel(obj.t));
            end

            % Validate BCs
            for s = {obj.bc, obj.bcx, obj.bcy, obj.bcz}
                if ~ismember(s{1},{'obc','pbc'})
                    error('UQMaT:Params:badBC','BC must be ''obc'' or ''pbc''.');
                end
            end

            % Validate disorder
            valid_dis = {'none','anderson','gaussian','binary', ...
                         'hop_disorder','impurity','vacancy'};
            if ~ismember(obj.disorder, valid_dis)
                error('UQMaT:Params:badDisorder', ...
                    'disorder must be: %s', strjoin(valid_dis,'|'));
            end

            % Auto-enable Anderson when W > 0
            if obj.W > 0 && strcmp(obj.disorder,'none')
                obj.disorder = 'anderson';
            end

            % Auto-enable sc from any pairing amplitude
            if any(obj.sc_Delta_s ~= 0) || obj.sc_Delta_t ~= 0 || any(obj.Delta_p ~= 0)
                obj.sc = true;
            end
            % pairing_mode
            has_s = any(obj.sc_Delta_s ~= 0);
            has_t = obj.sc_Delta_t ~= 0;
            if has_s && has_t,  obj.pairing_mode = 'mixed';
            elseif has_t,       obj.pairing_mode = 'triplet';
            else,               obj.pairing_mode = 'singlet';
            end

            % Spin flag — also enabled by altermagnet (needs σ_z)
            has_B  = any(obj.B ~= 0);
            has_R  = obj.lambda_R ~= 0;
            has_AM = ~isempty(obj.am_components);
            obj.spin = obj.sc || has_B || has_R || has_AM;

            % On-site energy
            obj.onsite_val = obj.onsite;

            % Expand Delta_s from SC struct
            N  = obj.N;
            Ns = N;
            if strcmp(obj.system,'kagome')
                Ns = 3 * obj.Nx * obj.Ny;
            elseif strcmp(obj.system,'tmd')
                Ns = 3 * obj.Nx * obj.Ny;
            elseif obj.is_2d
                Ns = obj.Nx * obj.Ny;
            elseif strcmp(obj.system,'diamond')
                Ns = 2 * obj.Nx * obj.Ny * obj.Nz;
            elseif obj.is_3d
                Ns = obj.Nx * obj.Ny * obj.Nz;
            end

            % Expand sc_Delta_s to Ns×1 site-resolved vector
            if isscalar(obj.sc_Delta_s)
                obj.Delta_s_vec = repmat(obj.sc_Delta_s, Ns, 1);
            elseif isvector(obj.sc_Delta_s) && numel(obj.sc_Delta_s) == Ns
                obj.Delta_s_vec = obj.sc_Delta_s(:);
            elseif isvector(obj.sc_Delta_s) && numel(obj.sc_Delta_s) ~= Ns
                error('UQMaT:Params:badDeltaS', ...
                    'SC.Delta_s vector length %d does not match Ns=%d.', ...
                    numel(obj.sc_Delta_s), Ns);
            else
                error('UQMaT:Params:badDeltaS', ...
                    'SC.Delta_s must be a scalar or Ns=%d vector.', Ns);
            end
            obj.Delta_t_val = obj.sc_Delta_t;

            % Expand Delta_p (Kitaev)
            Nb = N - 1;
            if isscalar(obj.Delta_p)
                obj.Delta_p_vec = repmat(obj.Delta_p, Nb, 1);
            elseif isvector(obj.Delta_p) && numel(obj.Delta_p) == Nb
                obj.Delta_p_vec = obj.Delta_p(:);
            else
                error('UQMaT:Params:badDeltaP', ...
                    'Delta_p: scalar or vector(len=N-1=%d). Got length %d.', Nb, numel(obj.Delta_p));
            end

            % Resolve impurity strengths
            if ~isempty(obj.imp_pos)
                if obj.is_2d
                    inp = obj.imp_pos(:);
                    Ni  = max(1, floor(numel(inp)/2));
                else
                    Ni = numel(obj.imp_pos);
                end
                if isscalar(obj.imp_potential)
                    obj.imp_potential_val = repmat(obj.imp_potential, Ni, 1);
                else
                    obj.imp_potential_val = obj.imp_potential(:);
                end
                if isscalar(obj.imp_magnetic)
                    obj.imp_magnetic_val = repmat(obj.imp_magnetic, Ni, 1);
                else
                    obj.imp_magnetic_val = obj.imp_magnetic(:);
                end
                if strcmp(obj.disorder,'none')
                    obj.disorder = 'impurity';
                end
            else
                obj.imp_potential_val = obj.imp_potential(:);
                obj.imp_magnetic_val  = obj.imp_magnetic(:);
            end
        end

        function disp(obj)
            fprintf('\n  UQMaT.Params\n  ============\n');
            fprintf('  system   : %s\n', obj.system);
            if obj.is_2d
                fprintf('  geometry : %d×%d = %d sites  bcx=%s  bcy=%s\n', ...
                    obj.Nx, obj.Ny, obj.Nx*obj.Ny, upper(obj.bcx), upper(obj.bcy));
            else
                fprintf('  geometry : N=%d  bc=%s\n', obj.N, upper(obj.bc));
            end
            if obj.tx == obj.ty
                fprintf('  hopping  : t=%.4g  mu=%.4g  t2=%.4g\n', obj.tx, obj.mu, obj.t2);
            else
                fprintf('  hopping  : tx=%.4g  ty=%.4g  mu=%.4g  t2=%.4g\n', ...
                    obj.tx, obj.ty, obj.mu, obj.t2);
            end
            if obj.onsite ~= 0
                fprintf('  onsite   : %.4g\n', obj.onsite_val);
            end
            if strcmp(obj.system,'ssh')
                fprintf('  SSH      : v=%.4g (intra)  w=%.4g (inter)\n', obj.v, obj.w);
            end
            if strcmp(obj.system,'kitaev')
                dp = mean(obj.Delta_p_vec);
                if dp ~= 0
                    fprintf('  Kitaev   : Delta_p=%.4g  →  spinless BdG 2N×2N\n', dp);
                else
                    fprintf('  Kitaev   : Delta_p=0  →  normal N×N\n');
                end
            end
            if any(obj.B ~= 0)
                fprintf('  Zeeman   : B=[%.4g, %.4g, %.4g]  (Bx,By,Bz)\n', ...
                    obj.B(1), obj.B(2), obj.B(3));
            end
            if obj.lambda_R ~= 0
                fprintf('  Rashba   : lambda_R=%.4g\n', obj.lambda_R);
            end
            if obj.t_nh ~= 0
                fprintf('  NH asym  : t_nh=%.4g  →  t_R=%.4g  t_L=%.4g\n', ...
                    obj.t_nh, obj.tx - obj.t_nh, obj.tx + obj.t_nh);
            end
            if obj.sc
                fprintf('  SC (BdG) : ON  Psi=(c_up,c_dn,c_up†,c_dn†)  4N×4N\n');
                fprintf('  pairing  : %s\n', obj.pairing_mode);
                if any(obj.sc_Delta_s ~= 0)
                    if isscalar(obj.sc_Delta_s)
                        fprintf('    Delta_s = %.4g  sym = %s\n', obj.sc_Delta_s, obj.sc_sym_s);
                    else
                        fprintf('    Delta_s = [vector, mean=%.4g, max=%.4g]  sym = %s\n', ...
                            mean(obj.sc_Delta_s), max(obj.sc_Delta_s), obj.sc_sym_s);
                    end
                end
                if obj.sc_Delta_t ~= 0
                    fprintf('    Delta_t = %.4g  sym_t = %s\n', obj.sc_Delta_t, obj.sc_sym_t);
                end
            else
                if obj.spin
                    fprintf('  SC       : off  (spinful 2N×2N)\n');
                else
                    fprintf('  SC       : off  (normal  N×N)\n');
                end
            end
            % Altermagnet
            for k = 1:numel(obj.am_components)
                c = obj.am_components{k};
                fprintf('  altermagnet: J=%.4g  sym=%s\n', c.J, c.sym);
            end
            % CDW
            if obj.cdw_Delta ~= 0
                fprintf('  CDW        : Delta=%.4g  sym=%s\n', obj.cdw_Delta, obj.cdw_sym);
            end
            switch obj.disorder
                case 'none'
                    fprintf('  disorder : none\n');
                case 'anderson'
                    fprintf('  disorder : Anderson  W=%.4g  seed=%d  (zero-mean uniform)\n', ...
                        obj.W, obj.seed);
                case 'gaussian'
                    fprintf('  disorder : Gaussian  W=%.4g  seed=%d  (zero-mean)\n', ...
                        obj.W, obj.seed);
                case 'binary'
                    fprintf('  disorder : binary alloy  W=%.4g  seed=%d\n', ...
                        obj.W, obj.seed);
                case 'hop_disorder'
                    fprintf('  disorder : hopping disorder  W=%.4g  seed=%d\n', ...
                        obj.W, obj.seed);
                case {'impurity','vacancy'}
                    if obj.is_2d
                        Ni = max(1, floor(numel(obj.imp_pos(:))/2));
                    else
                        Ni = numel(obj.imp_pos);
                    end
                    fprintf('  disorder : %s  Ni=%d  sigma=%.3g  geometry=%s\n', ...
                        obj.disorder, Ni, obj.imp_sigma, obj.imp_geometry);
                    if any(obj.imp_potential_val ~= 0)
                        fprintf('    V0 = %s\n', num2str(obj.imp_potential_val.','%.4g '));
                    end
                    if any(obj.imp_magnetic_val ~= 0)
                        fprintf('    J  = %s\n', num2str(obj.imp_magnetic_val.','%.4g '));
                    end
            end
            fprintf('  broadening : eta=%.4g\n', obj.eta);
            if ~isempty(obj.EnergyGrid)
                fprintf('  EnergyGrid : %d points  [%.4g, %.4g]  step≈%.4g\n\n', ...
                    numel(obj.EnergyGrid), obj.EnergyGrid(1), obj.EnergyGrid(end), ...
                    mean(diff(obj.EnergyGrid)));
            else
                fprintf('  EnergyGrid : auto (from spectrum)\n\n');
            end
        end
    end

    methods (Access = private)

        function obj = parseSC(obj)
            % Parse SC struct into resolved sc_Delta_s, sc_sym_s, sc_Delta_t, sc_sym_t
            % SC = struct('Delta_s',0.3,'sym','onsite','Delta_t',0.1,'sym_t','px')
            sc = obj.SC;
            if isempty(fieldnames(sc)), return; end

            if isfield(sc,'Delta_s'), obj.sc_Delta_s = sc.Delta_s; end
            if isfield(sc,'sym'),     obj.sc_sym_s   = lower(sc.sym); end
            if isfield(sc,'sym_s'),   obj.sc_sym_s   = lower(sc.sym_s); end
            if isfield(sc,'Delta_t'), obj.sc_Delta_t = sc.Delta_t; end
            if isfield(sc,'sym_t'),   obj.sc_sym_t   = lower(sc.sym_t); end
        end

        function obj = parseAM(obj)
            % Parse AM struct into am_components cell array.
            % Simple form:
            %   AM = struct('J',0.5,'sym','gwave')
            % Multiple components — any field named J_<label> with sym_<label>:
            %   AM = struct('J_s',0.3,'sym_s','onsite','J_dxy',0.2,'sym_dxy','dxy')
            % Or a cell array of structs:
            %   AM = {struct('J',0.3,'sym','onsite'), struct('J',0.2,'sym','dxy')}
            am = obj.AM;
            if isempty(am), return; end

            comps = {};

            if iscell(am)
                % Cell array of structs
                for k = 1:numel(am)
                    c = am{k};
                    comps{end+1} = struct('J', c.J, 'sym', lower(c.sym)); %#ok<AGROW>
                end

            elseif isstruct(am)
                f = fieldnames(am);

                % Check for simple J + sym form
                if isfield(am,'J') && isfield(am,'sym')
                    comps{end+1} = struct('J', am.J, 'sym', lower(am.sym));
                end

                % Check for J_<label> + sym_<label> pairs
                % e.g. J_dxy → look for sym_dxy
                for k = 1:numel(f)
                    fname = f{k};
                    if length(fname) > 2 && strncmp(fname,'J_',2)
                        label = fname(3:end);   % e.g. 'dxy' from 'J_dxy'
                        sym_field = ['sym_' label];
                        J_val = am.(fname);
                        if isfield(am, sym_field)
                            sym_val = lower(am.(sym_field));
                        else
                            sym_val = label;   % use the label itself as sym
                        end
                        comps{end+1} = struct('J', J_val, 'sym', sym_val); %#ok<AGROW>
                    end
                end
            end

            obj.am_components = comps;
        end

        function obj = parseCDW(obj)
            % Parse CDW struct into cdw_Delta and cdw_sym.
            % CDW = struct('Delta',0.2,'sym','onsite')
            cdw = obj.CDW;
            if isempty(fieldnames(cdw)), return; end

            if isfield(cdw,'Delta'), obj.cdw_Delta = cdw.Delta; end
            if isfield(cdw,'sym'),   obj.cdw_sym   = lower(cdw.sym); end
        end

    end

    methods (Static, Access = private)
        function s = paramList()
            s = ['  system, N, bc, Nx, Ny, bcx, bcy,\n' ...
                 '  t, mu, onsite, t2, v, w,\n' ...
                 '  B, lambda_R, t_nh (non-Hermitian asymmetry),\n' ...
                 '  SC (struct: Delta_s,sym,Delta_t,sym_t),\n' ...
                 '  AM (struct: J,sym or J_<sym>,sym_<sym>,...),\n' ...
                 '  CDW (struct: Delta,sym), Delta_p,\n' ...
                 '  disorder, W, seed,\n' ...
                 '  imp_pos, imp_potential, imp_magnetic, imp_sigma, imp_geometry,\n' ...
                 '  eta, EnergyGrid'];
        end
    end
end
