classdef Disorder
    % UQMaT.Observables.Conventional.Disorder
    %   Disorder-averaged spectral quantities via ensemble sampling.
    %
    %   For a single realisation, just use UQMaT.Observables.Conventional.Spectral.
    %   This class runs Nreal disorder realisations and averages DOS / LDOS.
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params('disorder','anderson','W',1.5,'N',80);
    %     obs = UQMaT.Observables.Conventional.Disorder(p)
    %     obs = UQMaT.Observables.Conventional.Disorder(p, 'Nreal', 200)
    %
    %     [w, D]      = obs.dos()           % disorder-averaged DOS
    %     [w, L]      = obs.ldos()          % averaged LDOS (Nomega × N)
    %     [w, L]      = obs.ldos(sites)
    %     ipr         = obs.ipr()           % inverse participation ratio per state
    %     [w, D, Dstd]= obs.dos()           % std dev as third output
    %
    %   Notes
    %   -----
    %   - Each realisation uses a different RNG seed (seed + realisation index).
    %   - IPR = sum_i |psi_i|^4 ; IPR → 1/N for extended, → 1 for localised.
    %   - The solver used per realisation is Direct (exact diag).

    properties
        Nreal (1,1) double {mustBeInteger,mustBePositive} = 100
    end

    properties (SetAccess = private)
        base_params UQMaT.Params
    end

    methods
        function obj = Disorder(p, varargin)
            if ~isa(p,'UQMaT.Params')
                error('UQMaT:Disorder:badInput', ...
                    'First argument must be a UQMaT.Params object.');
            end
            if strcmp(p.disorder,'none')
                warning('UQMaT:Disorder:noDisorder', ...
                    'disorder=''none'' — ensemble averaging has no effect.');
            end
            obj.base_params = p;
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        % ---- Disorder-averaged DOS -------------------------------------
        function [omega, D, Dstd] = dos(obj, omega_in)
            p0  = obj.base_params;
            sys0 = UQMaT.System(p0);
            E0   = sys0.E;
            if nargin < 2 || isempty(omega_in)
                pad   = max(0.2*(E0(end)-E0(1)), 4*p0.eta);
                omega = linspace(E0(1)-pad, E0(end)+pad, p0.Nomega)';
            else
                omega = omega_in(:);
            end

            Nw   = numel(omega);
            Dall = zeros(Nw, obj.Nreal);

            for ir = 1:obj.Nreal
                p_r  = obj.seedParams(p0, ir);
                sys  = UQMaT.System(p_r);
                slv  = UQMaT.Solver.Diagonalization.Direct(sys);
                obs  = UQMaT.Observables.Conventional.Spectral(slv);
                [~, Dr] = obs.dos('omega', omega);
                Dall(:,ir) = Dr;
            end

            D    = mean(Dall, 2);
            Dstd = std(Dall, 0, 2);
        end

        % ---- Disorder-averaged LDOS ------------------------------------
        function [omega, L] = ldos(obj, sites, omega_in)
            p0  = obj.base_params;
            Nph = UQMaT.System(p0).Ndim;
            if nargin < 2 || isempty(sites), sites = 1:Nph; end
            if nargin < 3, omega_in = []; end

            sys0 = UQMaT.System(p0);
            E0   = sys0.E;
            if isempty(omega_in)
                pad   = max(0.2*(E0(end)-E0(1)), 4*p0.eta);
                omega = linspace(E0(1)-pad, E0(end)+pad, p0.Nomega)';
            else
                omega = omega_in(:);
            end

            Nw   = numel(omega);
            Ns   = numel(sites);
            Lsum = zeros(Nw, Ns);

            for ir = 1:obj.Nreal
                p_r  = obj.seedParams(p0, ir);
                sys  = UQMaT.System(p_r);
                slv  = UQMaT.Solver.Diagonalization.Direct(sys);
                obs  = UQMaT.Observables.Conventional.Spectral(slv);
                [~, Lr] = obs.ldos(sites, 'omega', omega);
                Lsum = Lsum + Lr;
            end

            L = Lsum / obj.Nreal;
        end

        % ---- Inverse Participation Ratio -------------------------------
        function [ipr, E] = ipr(obj)
            % IPR = sum_i |psi_i|^4  averaged over disorder realisations.
            % Returns ipr (Ndim_H × 1) and corresponding eigenvalues E.
            p0   = obj.base_params;
            sys0 = UQMaT.System(p0);
            Nst  = numel(sys0.E);
            Nph  = sys0.Ndim;

            ipr_sum = zeros(Nst, 1);
            E_sum   = zeros(Nst, 1);

            for ir = 1:obj.Nreal
                p_r = obj.seedParams(p0, ir);
                sys = UQMaT.System(p_r);
                V   = sys.V(1:Nph, :);   % electron block
                ipr_r = sum(abs(V).^4, 1)';
                ipr_sum = ipr_sum + ipr_r;
                E_sum   = E_sum + sys.E;
            end

            ipr = ipr_sum / obj.Nreal;
            E   = E_sum   / obj.Nreal;
        end

        function disp(obj)
            fprintf('  UQMaT.Observables.Conventional.Disorder\n');
            fprintf('  Model : %s, W=%.4g, Nreal=%d\n', ...
                upper(obj.base_params.model), ...
                obj.base_params.W, obj.Nreal);
        end
    end

    methods (Access = private)
        function p_r = seedParams(~, p0, ir)
            % Return a copy of p0 with a unique seed per realisation
            p_r = UQMaT.Params( ...
                'model',    p0.model, ...
                'N',        p0.N, ...
                'bc',       p0.bc, ...
                't',        p0.t, ...
                'mu',       p0.mu, ...
                't2',       p0.t2, ...
                'v',        p0.v, ...
                'w',        p0.w, ...
                'bdg',      p0.bdg, ...
                'disorder', p0.disorder, ...
                'W',        p0.W, ...
                'seed',     p0.seed + ir, ...
                'eta',      p0.eta, ...
                'Nomega',   p0.Nomega);
        end
    end
end
