classdef ParamsKKRGreen
    % UQMaT.ParamsKKRGreen
    %   Parameters for the KKR (Korringa-Kohn-Rostoker) Green's function solver.
    %
    %   Computes G via the KKR secular equation:
    %     G^{KKR}(ω) = [t(ω)^{-1} - B]^{-1}
    %   where t = single-site t-matrix, B = structure constants.
    %
    %   Usage:
    %     p   = UQMaT.ParamsKKRGreen('system','square','Nx',10,'Ny',10, ...
    %               'EnergyGrid',-3:0.01:3, 'eta',0.01)
    %     sys = UQMaT.System(p.toParams())
    %     slv = UQMaT.Solver.Greens.KKR(sys, p)
    %     ev  = UQMaT.EV(slv)
    %     [w, D]   = ev.dos()
    %     [w, L]   = slv.ldos(sites)
    %     [w, Sig] = slv.self_energy()   % CPA self-energy (use_cpa=true)
    %
    %   Key parameters:
    %     use_cpa  : enable Coherent Potential Approximation for disorder
    %     cpa_tol  : CPA convergence threshold
    %
    %  Citation
    %  --------
    %  Korringa J., Physica 13, 392 (1947).
    %  Kohn W., Rostoker N., Phys. Rev. 94, 1111 (1954).
    %  Faulkner J.S., Stocks G.M., Phys. Rev. B 21, 3222 (1980).

    properties
        system    (1,:) char   = 'chain'
        N         (1,1) double {mustBePositive,mustBeInteger} = 50
        Nx        (1,1) double {mustBePositive,mustBeInteger} = 10
        Ny        (1,1) double {mustBePositive,mustBeInteger} = 10
        bc        (1,:) char   = 'obc'
        bcx       (1,:) char   = 'obc'
        bcy       (1,:) char   = 'obc'
        t         double       = 1.0
        mu        (1,1) double = 0.0
        SC        struct       = struct()
        W         (1,1) double = 0.0
        disorder  (1,:) char   = 'none'
        seed      (1,1) double = 0

        % KKR-specific
        EnergyGrid (1,:) double = []
        eta        (1,1) double {mustBePositive} = 0.01
        use_cpa    (1,1) logical = false
        cpa_tol    (1,1) double  = 1e-6
        cpa_maxit  (1,1) double  = 100
    end

    methods
        function obj = ParamsKKRGreen(varargin)
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
                't',obj.t,'mu',obj.mu,'SC',obj.SC, ...
                'W',obj.W,'disorder',obj.disorder,'seed',obj.seed, ...
                'eta',obj.eta,'EnergyGrid',obj.EnergyGrid);
        end

        function [ev, slv] = solve(obj)
            sys = UQMaT.System(obj.toParams());
            slv = UQMaT.Solver.Greens.KKR(sys, obj);
            ev  = UQMaT.EV(slv);
        end

        function disp(obj)
            fprintf('\n  UQMaT.ParamsKKRGreen\n');
            fprintf('  system=%s  N=%d  eta=%.4g  use_cpa=%d\n', ...
                obj.system, obj.N, obj.eta, obj.use_cpa);
            fprintf('  EnergyGrid: %d points [%.4g, %.4g]\n\n', ...
                numel(obj.EnergyGrid), obj.EnergyGrid(1), obj.EnergyGrid(end));
        end
    end
end
