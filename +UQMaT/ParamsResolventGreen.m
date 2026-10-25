classdef ParamsResolventGreen
    % UQMaT.ParamsResolventGreen
    %   Parameters for the resolvent (direct inversion) Green's function solver.
    %
    %   G^r(ω) = (ω + iη - H)^{-1}   computed by direct matrix inversion.
    %   Exact but O(N³) per energy point — best for moderate N (≲ 1000).
    %
    %   Usage:
    %     p   = UQMaT.ParamsResolventGreen('system','chain','N',200, ...
    %               'EnergyGrid',-3:0.005:3, 'eta',0.01)
    %     sys = UQMaT.System(p.toParams())
    %     slv = UQMaT.Solver.Greens.Resolvent(sys)
    %     ev  = UQMaT.EV(slv)
    %     [w, D] = ev.dos()
    %
    %   Key parameters:
    %     EnergyGrid : energy points ω (required — no auto fallback for GF solvers)
    %     eta        : broadening η = Im(ω + iη) — controls peak width and stability
    %                  Smaller η → sharper features; too small → numerical noise
    %     use_sparse : use sparse matrix inversion (faster for large sparse H)
    %
    %  Citation
    %  --------
    %  Economou E.N., Green's Functions in Quantum Physics, Springer (2006).

    properties
        % System (passed to UQMaT.Params)
        system    (1,:) char   = 'chain'
        N         (1,1) double {mustBePositive,mustBeInteger} = 50
        Nx        (1,1) double {mustBePositive,mustBeInteger} = 10
        Ny        (1,1) double {mustBePositive,mustBeInteger} = 10
        Nz        (1,1) double {mustBePositive,mustBeInteger} = 10
        bc        (1,:) char   = 'obc'
        bcx       (1,:) char   = 'obc'
        bcy       (1,:) char   = 'obc'
        t         double       = 1.0
        mu        (1,1) double = 0.0
        SC        struct       = struct()
        AM        struct       = struct()
        CDW       struct       = struct()
        Delta_p   double       = 0.0
        B         (1,3) double = [0,0,0]
        lambda_R  (1,1) double = 0.0
        W         (1,1) double = 0.0
        disorder  (1,:) char   = 'none'
        seed      (1,1) double = 0

        % Resolvent-specific
        EnergyGrid (1,:) double = []   % REQUIRED: energy grid
        eta        (1,1) double {mustBePositive} = 0.01
        use_sparse (1,1) logical = false
    end

    methods
        function obj = ParamsResolventGreen(varargin)
            for k = 1:2:numel(varargin)
                obj.(varargin{k}) = varargin{k+1};
            end
            if isempty(obj.EnergyGrid)
                warning('UQMaT:PResolvent:noGrid', ...
                    'EnergyGrid not set. Using auto range [-4,4] with 500 points.');
                obj.EnergyGrid = linspace(-4, 4, 500);
            end
        end

        function p = toParams(obj)
            p = UQMaT.Params('system',obj.system,'N',obj.N, ...
                'Nx',obj.Nx,'Ny',obj.Ny,'Nz',obj.Nz, ...
                'bc',obj.bc,'bcx',obj.bcx,'bcy',obj.bcy, ...
                't',obj.t,'mu',obj.mu,'SC',obj.SC,'AM',obj.AM,'CDW',obj.CDW, ...
                'Delta_p',obj.Delta_p,'B',obj.B,'lambda_R',obj.lambda_R, ...
                'W',obj.W,'disorder',obj.disorder,'seed',obj.seed, ...
                'eta',obj.eta,'EnergyGrid',obj.EnergyGrid);
        end

        function ev = solve(obj)
            sys = UQMaT.System(obj.toParams());
            ev  = UQMaT.EV(sys, 'resolvent', 'eta', obj.eta);
        end

        function disp(obj)
            fprintf('\n  UQMaT.ParamsResolventGreen\n');
            fprintf('  system=%s  N=%d  eta=%.4g\n', obj.system, obj.N, obj.eta);
            fprintf('  EnergyGrid: %d points [%.4g, %.4g]\n', ...
                numel(obj.EnergyGrid), obj.EnergyGrid(1), obj.EnergyGrid(end));
            fprintf('  use_sparse=%d\n\n', obj.use_sparse);
        end
    end
end
