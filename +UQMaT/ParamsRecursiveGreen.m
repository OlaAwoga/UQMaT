classdef ParamsRecursiveGreen
    % UQMaT.ParamsRecursiveGreen
    %   Parameters for the recursive Green's function (RGF) solver.
    %
    %   Propagates the surface GF layer by layer:
    %     G^r(n) = [ω - H_n - Σ_L(n)]^{-1}
    %   Scales as O(N · M³) where M = layer block size.
    %   Memory O(M²) — ideal for long quasi-1D systems and transport.
    %
    %   Usage:
    %     p   = UQMaT.ParamsRecursiveGreen('system','chain','N',1000, ...
    %               'EnergyGrid',-3:0.01:3, 'eta',0.01, 'n_layers',1000)
    %     sys = UQMaT.System(p.toParams())
    %     slv = UQMaT.Solver.Greens.Recursive(sys, p)
    %     ev  = UQMaT.EV(slv)
    %     [w, D] = ev.dos()
    %     [w, T] = slv.transmission()
    %
    %   Key parameters:
    %     n_layers : number of principal layers (0 = auto: sqrt(Hdim))
    %                More layers → smaller M → faster per layer, more layers
    %     eta      : broadening (also acts as lead coupling in WBL)
    %
    %  Citation
    %  --------
    %  López Sancho M.P. et al., J. Phys. F 14, 1205 (1984).
    %  Ferry D.K., Goodnick S.M., Transport in Nanostructures (1997).

    properties
        system    (1,:) char   = 'chain'
        N         (1,1) double {mustBePositive,mustBeInteger} = 200
        Nx        (1,1) double {mustBePositive,mustBeInteger} = 10
        Ny        (1,1) double {mustBePositive,mustBeInteger} = 10
        bc        (1,:) char   = 'obc'
        bcx       (1,:) char   = 'obc'
        bcy       (1,:) char   = 'obc'
        t         double       = 1.0
        mu        (1,1) double = 0.0
        SC        struct       = struct()
        AM        struct       = struct()
        B         (1,3) double = [0,0,0]
        lambda_R  (1,1) double = 0.0
        W         (1,1) double = 0.0
        disorder  (1,:) char   = 'none'
        seed      (1,1) double = 0

        % RGF-specific
        EnergyGrid (1,:) double = []
        eta        (1,1) double {mustBePositive} = 0.01
        n_layers   (1,1) double {mustBeInteger} = 0   % 0 = auto, >0 = user-specified
    end

    methods
        function obj = ParamsRecursiveGreen(varargin)
            for k = 1:2:numel(varargin)
                obj.(varargin{k}) = varargin{k+1};
            end
            if isempty(obj.EnergyGrid)
                obj.EnergyGrid = linspace(-4, 4, 500);
            end
        end

        function p = toParams(obj)
            p = UQMaT.Params('system',obj.system,'N',obj.N, ...
                'Nx',obj.Nx,'Ny',obj.Ny,'bc',obj.bc,'bcx',obj.bcx,'bcy',obj.bcy, ...
                't',obj.t,'mu',obj.mu,'SC',obj.SC,'AM',obj.AM, ...
                'B',obj.B,'lambda_R',obj.lambda_R, ...
                'W',obj.W,'disorder',obj.disorder,'seed',obj.seed, ...
                'eta',obj.eta,'EnergyGrid',obj.EnergyGrid);
        end

        function [ev, slv] = solve(obj)
            sys = UQMaT.System(obj.toParams());
            slv = UQMaT.Solver.Greens.Recursive(sys, obj);
            ev  = UQMaT.EV(slv);
        end

        function disp(obj)
            fprintf('\n  UQMaT.ParamsRecursiveGreen\n');
            fprintf('  system=%s  N=%d  n_layers=%d  eta=%.4g\n', ...
                obj.system, obj.N, obj.n_layers, obj.eta);
            fprintf('  EnergyGrid: %d points [%.4g, %.4g]\n\n', ...
                numel(obj.EnergyGrid), obj.EnergyGrid(1), obj.EnergyGrid(end));
        end
    end
end
