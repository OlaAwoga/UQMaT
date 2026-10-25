classdef ParamsChebyshev
    % UQMaT.ParamsChebyshev
    %   Parameters for the Chebyshev / Kernel Polynomial Method (KPM) solver.
    %
    %   Expands DOS as a Chebyshev series:
    %     D(ω) = 1/(π√(1-ω²)) · Σ_m μ_m T_m(ω)
    %   where μ_m = Tr[T_m(H)] are moments computed via matrix-vector products.
    %   Scales as O(N·M·Nrand) — no matrix inversion, ideal for very large N.
    %
    %   Usage:
    %     p   = UQMaT.ParamsChebyshev('system','square','Nx',50,'Ny',50, ...
    %               'M',1024, 'Nrand',20, 'EnergyGrid',-3:0.005:3)
    %     sys = UQMaT.System(p.toParams())
    %     ev  = UQMaT.EV(sys, 'chebyshev', 'M',p.M, 'Nrand',p.Nrand)
    %     [w, D] = ev.dos()
    %
    %   Key parameters:
    %     M      : number of Chebyshev moments (resolution ~ π/M in scaled units)
    %              More moments → sharper features. Typical: 256-4096.
    %     Nrand  : number of stochastic trace vectors (more → less noise)
    %     kernel : 'jackson' (smooth, default) | 'lorentz' (Lorentzian broadening)
    %     eta    : broadening for Lorentz kernel (ignored for Jackson kernel)
    %
    %  Citation
    %  --------
    %  Weisse A. et al., Rev. Mod. Phys. 78, 275 (2006). [KPM review]
    %  Silver R.N., Roder H., Phys. Rev. E 56, 4822 (1997). [stochastic trace]

    properties
        system    (1,:) char   = 'square'
        N         (1,1) double {mustBePositive,mustBeInteger} = 50
        Nx        (1,1) double {mustBePositive,mustBeInteger} = 30
        Ny        (1,1) double {mustBePositive,mustBeInteger} = 30
        Nz        (1,1) double {mustBePositive,mustBeInteger} = 10
        bc        (1,:) char   = 'obc'
        bcx       (1,:) char   = 'obc'
        bcy       (1,:) char   = 'obc'
        t         double       = 1.0
        mu        (1,1) double = 0.0
        SC        struct       = struct()
        AM        struct       = struct()
        W         (1,1) double = 0.0
        disorder  (1,:) char   = 'none'
        seed      (1,1) double = 0

        % KPM-specific
        EnergyGrid (1,:) double = []
        M          (1,1) double {mustBePositive,mustBeInteger} = 512
        Nrand      (1,1) double {mustBePositive,mustBeInteger} = 20
        kernel     (1,:) char   = 'jackson'
        eta        (1,1) double {mustBePositive} = 0.01
        scale_factor (1,1) double = 1.05  % H rescaled to [-1,1]/scale_factor
    end

    methods
        function obj = ParamsChebyshev(varargin)
            for k = 1:2:numel(varargin)
                obj.(varargin{k}) = varargin{k+1};
            end
            if isempty(obj.EnergyGrid)
                obj.EnergyGrid = linspace(-4, 4, 600);
            end
        end

        function p = toParams(obj)
            p = UQMaT.Params('system',obj.system,'N',obj.N, ...
                'Nx',obj.Nx,'Ny',obj.Ny,'Nz',obj.Nz, ...
                'bc',obj.bc,'bcx',obj.bcx,'bcy',obj.bcy, ...
                't',obj.t,'mu',obj.mu,'SC',obj.SC,'AM',obj.AM, ...
                'W',obj.W,'disorder',obj.disorder,'seed',obj.seed, ...
                'eta',obj.eta,'EnergyGrid',obj.EnergyGrid);
        end

        function ev = solve(obj)
            sys = UQMaT.System(obj.toParams());
            ev  = UQMaT.EV(sys, 'chebyshev', ...
                'M', obj.M, 'Nrand', obj.Nrand, 'kernel', obj.kernel);
        end

        function disp(obj)
            fprintf('\n  UQMaT.ParamsChebyshev (KPM)\n');
            fprintf('  system=%s  Nx=%d  Ny=%d\n', obj.system, obj.Nx, obj.Ny);
            fprintf('  M=%d  Nrand=%d  kernel=%s  eta=%.4g\n', ...
                obj.M, obj.Nrand, obj.kernel, obj.eta);
            fprintf('  EnergyGrid: %d points [%.4g, %.4g]\n\n', ...
                numel(obj.EnergyGrid), obj.EnergyGrid(1), obj.EnergyGrid(end));
        end
    end
end
