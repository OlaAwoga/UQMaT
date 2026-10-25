classdef Params2D
    % UQMaT.Params2D  Parameters for 2D lattice systems.
    %
    %   A convenience wrapper for 2D systems that pre-sets sensible
    %   defaults for Nx, Ny, bcx, bcy and validates 2D-specific fields.
    %   Internally converted to UQMaT.Params for use with System and EV.
    %
    %   p = UQMaT.Params2D()
    %   p = UQMaT.Params2D('Nx',20,'Ny',20,'t',1,'mu',0)
    %   sys = UQMaT.System(p.toParams())
    %   ev  = UQMaT.EV(UQMaT.System(p.toParams()))
    %
    % ====================================================================
    %  GEOMETRY
    % ====================================================================
    %   Nx, Ny   : sites per direction              (default: 10, 10)
    %   bcx, bcy : BC per axis  'obc'|'pbc'        (default: 'obc','obc')
    %   system   : '2d lattice type'               (default: 'square')
    %              'square' | 'kagome' | 'tmd'
    %
    % ====================================================================
    %  NORMAL STATE
    % ====================================================================
    %   t        : hopping (scalar or [tx,ty])     (default: 1.0)
    %   mu       : chemical potential              (default: 0.0)
    %   t2       : next-nearest hopping           (default: 0.0)
    %   onsite   : uniform on-site energy         (default: 0.0)
    %
    % ====================================================================
    %  ZEEMAN / RASHBA
    % ====================================================================
    %   B        : Zeeman field [Bx,By,Bz]        (default: [0,0,0])
    %   lambda_R : Rashba SOC                     (default: 0.0)
    %
    % ====================================================================
    %  ORDER PARAMETERS (structs)
    % ====================================================================
    %   SC  = struct('Delta_s',0.3,'sym','onsite')        singlet SC (uniform)
    %   SC  = struct('Delta_s',Delta_vec,'sym','onsite')  site-resolved Δ(i) (Ns×1)
    %   SC  = struct('Delta_t',0.1,'sym_t','chiral_p')    triplet SC
    %   AM  = struct('J',0.5,'sym','gwave')          altermagnet
    %   CDW = struct('Delta',0.2,'sym','onsite')      CDW
    %   Delta_p : Kitaev p-wave (scalar)            (default: 0.0)
    %
    % ====================================================================
    %  DISORDER
    % ====================================================================
    %   W             : Anderson disorder          (default: 0.0)
    %   disorder      : type (auto 'anderson' when W>0)
    %   imp_pos       : impurity positions [ix,iy] or Nimp×2
    %   imp_potential : on-site potential          (default: 0.0)
    %   imp_magnetic  : exchange strength          (default: 0.0)
    %   imp_sigma     : Gaussian spread            (default: 0.0)
    %   imp_geometry  : 'point'|'hline'|'vline'|'line'|'cluster'
    %   seed          : RNG seed                   (default: 0)
    %
    % ====================================================================
    %  TMD
    % ====================================================================
    %   tmd_material : 'MoS2'|'WS2'|'MoSe2'|'WSe2'|'custom'
    %   tmd_params   : custom struct
    %
    % ====================================================================
    %  SPECTRAL
    % ====================================================================
    %   eta      : Lorentzian broadening           (default: 0.02)
    %   EnergyGrid : energy grid vector            (default: [])
    %
    %  Citation
    %  --------
    %  UQMaT: see CITATION.cff.

    properties
        % Geometry
        system   (1,:) char  = 'square'
        Nx       (1,1) double {mustBeInteger,mustBePositive} = 10
        Ny       (1,1) double {mustBeInteger,mustBePositive} = 10
        bcx      (1,:) char  = 'obc'
        bcy      (1,:) char  = 'obc'

        % Normal state
        t        double      = 1.0
        mu       (1,1) double = 0.0
        t2       (1,1) double = 0.0
        onsite   (1,1) double = 0.0

        % Zeeman / Rashba
        B        (1,3) double = [0 0 0]
        lambda_R (1,1) double = 0.0

        % Order parameters (structs — see Params)
        SC       struct = struct()
        AM       struct = struct()
        CDW      struct = struct()
        Delta_p  double = 0.0

        % Disorder
        W             (1,1) double = 0.0
        disorder      (1,:) char   = 'none'
        seed          (1,1) double {mustBeInteger,mustBeNonnegative} = 0
        imp_pos       double = []
        imp_potential double = 0.0
        imp_magnetic  double = 0.0
        imp_sigma     (1,1) double = 0.0
        imp_geometry  (1,:) char   = 'point'

        % TMD
        tmd_material (1,:) char = 'MoS2'
        tmd_params   struct = struct()

        % Spectral
        eta        (1,1) double {mustBePositive} = 0.02
        EnergyGrid (1,:) double = []
    end

    properties (SetAccess = private)
        Ns (1,1) double = 100   % Nx * Ny (or 3*Nx*Ny for kagome/tmd)
    end

    methods

        function obj = Params2D(varargin)
            if mod(numel(varargin),2) ~= 0
                error('UQMaT:Params2D:badInput','Arguments must be Name-Value pairs.');
            end
            for k = 1:2:numel(varargin)
                name = varargin{k};
                val  = varargin{k+1};
                if ~isprop(obj, name)
                    error('UQMaT:Params2D:unknownParam', ...
                        'Unknown parameter: "%s"\nSee UQMaT.Params2D for valid fields.', name);
                end
                obj.(name) = val;
            end

            % Normalise
            obj.system       = lower(strtrim(obj.system));
            obj.bcx          = lower(strtrim(obj.bcx));
            obj.bcy          = lower(strtrim(obj.bcy));
            obj.disorder     = lower(strtrim(obj.disorder));
            obj.imp_geometry = lower(strtrim(obj.imp_geometry));

            % Validate system
            valid_2d = {'square','kagome','tmd'};
            if ~ismember(obj.system, valid_2d)
                error('UQMaT:Params2D:badSystem', ...
                    'Params2D system must be: %s', strjoin(valid_2d,'|'));
            end

            % Validate BCs
            for bc_ = {obj.bcx, obj.bcy}
                if ~ismember(bc_{1},{'obc','pbc'})
                    error('UQMaT:Params2D:badBC','bcx/bcy must be ''obc'' or ''pbc''.');
                end
            end

            % Derived: Ns
            switch obj.system
                case 'square'
                    obj.Ns = obj.Nx * obj.Ny;
                case {'kagome','tmd'}
                    obj.Ns = 3 * obj.Nx * obj.Ny;
            end
        end

        function p = toParams(obj)
            % Convert to UQMaT.Params for use with System and EV.
            %   sys = UQMaT.System(p2d.toParams())
            %   ev  = UQMaT.EV(UQMaT.System(p2d.toParams()))
            p = UQMaT.Params(...
                'system',   obj.system, ...
                'Nx',       obj.Nx, ...
                'Ny',       obj.Ny, ...
                'bcx',      obj.bcx, ...
                'bcy',      obj.bcy, ...
                't',        obj.t, ...
                'mu',       obj.mu, ...
                't2',       obj.t2, ...
                'onsite',   obj.onsite, ...
                'B',        obj.B, ...
                'lambda_R', obj.lambda_R, ...
                'SC',       obj.SC, ...
                'AM',       obj.AM, ...
                'CDW',      obj.CDW, ...
                'Delta_p',  obj.Delta_p, ...
                'W',        obj.W, ...
                'disorder', obj.disorder, ...
                'seed',     obj.seed, ...
                'imp_pos',       obj.imp_pos, ...
                'imp_potential', obj.imp_potential, ...
                'imp_magnetic',  obj.imp_magnetic, ...
                'imp_sigma',     obj.imp_sigma, ...
                'imp_geometry',  obj.imp_geometry, ...
                'tmd_material',  obj.tmd_material, ...
                'tmd_params',    obj.tmd_params, ...
                'eta',           obj.eta, ...
                'EnergyGrid',    obj.EnergyGrid);
        end

        function ev = solve(obj)
            % Shortcut: build System and diagonalise in one call.
            %   ev = p2d.solve()
            %   [w, D] = ev.dos()
            ev = UQMaT.EV(UQMaT.System(obj.toParams()));
        end

        function disp(obj)
            fprintf('\n  UQMaT.Params2D\n');
            fprintf('  ==============\n');
            fprintf('  system   : %s\n', obj.system);
            fprintf('  geometry : %d × %d = %d sites   bcx=%s  bcy=%s\n', ...
                obj.Nx, obj.Ny, obj.Ns, upper(obj.bcx), upper(obj.bcy));
            fprintf('  t=%.4g  mu=%.4g  t2=%.4g\n', mean(obj.t), obj.mu, obj.t2);
            if ~isempty(fieldnames(obj.SC))
                fprintf('  SC       : ');
                f = fieldnames(obj.SC);
                for k = 1:numel(f)
                    val = obj.SC.(f{k});
                    if isnumeric(val) && ~isscalar(val)
                        fprintf('%s=[vector,mean=%.4g,max=%.4g] ', ...
                            f{k}, mean(val(:)), max(abs(val(:))));
                    elseif isnumeric(val)
                        fprintf('%s=%.4g ', f{k}, val);
                    else
                        fprintf('%s=%s ', f{k}, val);
                    end
                end
                fprintf('\n');
            end
            if ~isempty(fieldnames(obj.AM))
                fprintf('  AM       : J=%.4g\n', obj.AM.J);
            end
            if any(obj.B ~= 0)
                fprintf('  B        : [%.4g %.4g %.4g]\n', obj.B);
            end
            if obj.lambda_R ~= 0
                fprintf('  lambda_R : %.4g\n', obj.lambda_R);
            end
            if obj.W > 0
                fprintf('  disorder : %s  W=%.4g\n', obj.disorder, obj.W);
            end
            fprintf('  eta=%.4g\n\n', obj.eta);
        end

    end
end
