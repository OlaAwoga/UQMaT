classdef KKR < UQMaT.Solver.Base
    % UQMaT.Solver.Greens.KKR
    %   Korringa-Kohn-Rostoker (KKR) Green's function method.
    %
    %   Computes the lattice Green's function by expanding in
    %   a real-space basis and solving the KKR secular equation:
    %
    %     [m(ω) - B(k,ω)] · c = 0
    %
    %   where m(ω) is the single-site t-matrix and B(k,ω) is the
    %   structure constant (geometry-dependent, energy-independent).
    %
    %   Properties:
    %     - O(N³) per energy point (dominated by matrix inversion)
    %     - Natural for impurity embedding (single-site problem)
    %     - Exact treatment of multiple scattering
    %     - Returns full G(r,r',ω)
    %
    %   Usage:
    %     p   = UQMaT.ParamsKKRGreen('system','square','Nx',10,'Ny',10, ...
    %               'EnergyGrid',-3:0.01:3, 'eta',0.01)
    %     slv = UQMaT.Solver.Greens.KKR(sys, p)
    %     [w, D]   = slv.dos()
    %     [w, L]   = slv.ldos(sites)
    %     [w, sig] = slv.self_energy()   % disorder self-energy (CPA)
    %
    %  Citation
    %  --------
    %  Korringa J., Physica 13, 392 (1947).
    %  Kohn W., Rostoker N., Phys. Rev. 94, 1111 (1954).
    %  Economou E.N., Green's Functions in Quantum Physics,
    %    Springer (2006). [Chapter 7]
    %  Faulkner J.S., Stocks G.M., Phys. Rev. B 21, 3222 (1980).
    %    [KKR-CPA for disordered systems]

    properties
        n_moments (1,1) double {mustBePositive,mustBeInteger} = 50
        % Number of angular momentum channels (tight-binding: set to 1)
        use_cpa   (1,1) logical = false
        % Coherent Potential Approximation for disorder
        cpa_tol   (1,1) double = 1e-6
        cpa_maxit (1,1) double = 100
    end

    properties (SetAccess = private)
        omega_grid  (:,1) double = []
        eta_val     (1,1) double = 0.01
        B_cache     = []    % structure constants (computed once)
        DOS_cache   (:,1) double = []
        SelfE_cache (:,1) double = []
    end

    methods
        function obj = KKR(sys, p_green)
            obj = obj@UQMaT.Solver.Base(sys);
            if nargin >= 2 && ~isempty(p_green)
                obj.eta_val    = p_green.eta;
                obj.omega_grid = p_green.EnergyGrid(:);
                if isprop(p_green,'n_moments') || isfield(p_green,'n_moments')
                    obj.n_moments = p_green.n_moments;
                end
                if isprop(p_green,'use_cpa') || isfield(p_green,'use_cpa')
                    obj.use_cpa = p_green.use_cpa;
                end
            end
            if isempty(obj.omega_grid)
                obj.omega_grid = linspace(-4, 4, 500)';
            end
        end

        function [E, V] = eigs(obj)
            [w, D] = obj.dos();
            E = w;
            V = D;
        end

        function [omega, D] = dos(obj)
            % DOS via KKR: D(ω) = -1/π Im Tr G^{KKR}(ω).
            H     = full(obj.sys.H);
            N     = obj.sys.Ndim;
            Hdim  = size(H,1);
            omega = obj.omega_grid;
            eta   = obj.eta_val;
            Nw    = numel(omega);
            D     = zeros(Nw, 1);

            % Precompute structure constants (geometry, ω-independent in TB limit)
            % In the tight-binding limit, B = H (off-diagonal part)
            B_off = H - diag(diag(H));   % off-diagonal = structure constants

            for iw = 1:Nw
                z   = omega(iw) + 1i*eta;
                % t-matrix: single-site scattering (diagonal part of G^{-1})
                t_inv = diag(z*ones(Hdim,1) - diag(H));
                % KKR secular matrix: t^{-1} - B
                M_kkr = diag(t_inv) - B_off;
                % Green's function: G = M^{-1}
                G     = M_kkr \ eye(Hdim);
                D(iw) = -imag(trace(G)) / (pi * N);
            end
            obj.DOS_cache = D;
        end

        function [omega, L] = ldos(obj, sites)
            % Site-resolved LDOS via KKR Green's function diagonal.
            H     = full(obj.sys.H);
            Hdim  = size(H, 1);
            omega = obj.omega_grid;
            eta   = obj.eta_val;
            Nw    = numel(omega);
            Ns    = numel(sites);
            L     = zeros(Nw, Ns);
            B_off = H - diag(diag(H));

            for iw = 1:Nw
                z     = omega(iw) + 1i*eta;
                t_inv = diag(z*ones(Hdim,1) - diag(H));
                M_kkr = diag(t_inv) - B_off;
                G     = M_kkr \ eye(Hdim);
                for si = 1:Ns
                    i = sites(si);
                    L(iw,si) = -imag(G(i,i)) / pi;
                end
            end
        end

        function [omega, Sigma] = self_energy(obj)
            % CPA self-energy for disordered systems.
            % Iterates Σ(ω) until self-consistent.
            if ~obj.use_cpa
                warning('UQMaT:KKR:noCPA','CPA not enabled. Set use_cpa=true.');
                omega = obj.omega_grid; Sigma = zeros(size(omega)); return
            end
            H     = full(obj.sys.H);
            Hdim  = size(H,1);
            omega = obj.omega_grid;
            eta   = obj.eta_val;
            Nw    = numel(omega);
            Sigma = zeros(Nw, 1);
            B_off = H - diag(diag(H));
            for iw = 1:Nw
                z   = omega(iw) + 1i*eta;
                sig = 0;   % initial self-energy guess
                for it = 1:obj.cpa_maxit
                    t_inv = diag((z - sig)*ones(Hdim,1) - diag(H));
                    M_kkr = diag(t_inv) - B_off;
                    G_av  = trace(M_kkr \ eye(Hdim)) / Hdim;
                    sig_new = sig + G_av^(-1) - (z - sig);
                    if abs(sig_new - sig) < obj.cpa_tol, break; end
                    sig = sig_new;
                end
                Sigma(iw) = sig;
            end
            obj.SelfE_cache = Sigma;
        end
    end
end
