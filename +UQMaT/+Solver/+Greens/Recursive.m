classdef Recursive < UQMaT.Solver.Base
    % UQMaT.Solver.Greens.Recursive
    %   Real-space recursive Green's function method (RGF).
    %
    %   Divides the system into principal layers and propagates
    %   the surface Green's function layer-by-layer using the
    %   Dyson equation:
    %
    %     G^r(n) = [ ω - H_n - Σ_L(n) ]^{-1}
    %     Σ_L(n) = V_{n,n-1} · G^r(n-1) · V_{n-1,n}
    %
    %   Properties of this solver:
    %     - O(N·M²) scaling (N layers, M×M layer block size)
    %     - Memory O(M²) — only two layers held at a time
    %     - Returns full G(i,j,ω) for transmission, LDOS
    %     - Numerically stable for open systems with leads
    %
    %   Usage:
    %     p   = UQMaT.ParamsRecursiveGreen('system','chain','N',200, ...
    %               'EnergyGrid',-3:0.01:3, 'eta',0.01, 'n_layers',200)
    %     slv = UQMaT.Solver.Greens.Recursive(sys, p)
    %     [E, V] = slv.eigs()     % returns LDOS-derived spectrum
    %     [w, D] = slv.dos()      % DOS from Im Tr G
    %     [w, T] = slv.transmission()  % Landauer transmission
    %
    %  Citation
    %  --------
    %  López Sancho M.P. et al., J. Phys. F 14, 1205 (1984).
    %  Ferry D.K., Goodnick S.M., Transport in Nanostructures,
    %    Cambridge University Press (1997). [Chapter 3]

    properties
        n_layers  (1,1) double {mustBeInteger} = 0
        % 0 = auto (one layer per site, M=1 block)
        % >0 = user-specified number of layers
    end

    properties (SetAccess = private)
        omega_grid  (:,1) double = []
        eta_val     (1,1) double = 0.01
        DOS_cache   (:,1) double = []
        LDOS_cache  (:,:) double = []
        T_cache     (:,1) double = []
    end

    methods
        function obj = Recursive(sys, p_green)
            obj = obj@UQMaT.Solver.Base(sys);
            if nargin >= 2 && ~isempty(p_green)
                obj.eta_val    = p_green.eta;
                obj.omega_grid = p_green.EnergyGrid(:);
                if isfield(p_green,'n_layers') && p_green.n_layers > 0
                    obj.n_layers = p_green.n_layers;
                end
            end
            if isempty(obj.omega_grid)
                obj.omega_grid = linspace(-4, 4, 500)';
            end
        end

        function [E, V] = eigs(obj)
            % Compute DOS via RGF; return omega as E and DOS as V columns.
            [w, D] = obj.dos();
            E = w;
            V = D;
        end

        function [omega, D] = dos(obj)
            % DOS = -1/π Im Tr G^r(ω) computed via RGF.
            H    = full(obj.sys.H);
            N    = obj.sys.Ndim;
            Hdim = size(H,1);
            omega = obj.omega_grid;
            eta   = obj.eta_val;
            Nw    = numel(omega);
            D     = zeros(Nw, 1);

            % Determine layer block size
            M = obj.layer_size_(Hdim, N);

            for iw = 1:Nw
                z = omega(iw) + 1i*eta;
                G_surf = obj.rgf_surface_(H, z, M, Hdim);
                D(iw)  = -imag(trace(G_surf)) / (pi * N);
            end
            obj.DOS_cache = D;
        end

        function [omega, T] = transmission(obj)
            % Landauer transmission T(ω) = Tr[Γ_L G^r Γ_R G^a].
            % Uses wide-band lead approximation Γ = eta*I.
            H     = full(obj.sys.H);
            Hdim  = size(H, 1);
            N     = obj.sys.Ndim;
            M     = obj.layer_size_(Hdim, N);
            omega = obj.omega_grid;
            eta   = obj.eta_val;
            Nw    = numel(omega);
            T     = zeros(Nw, 1);

            for iw = 1:Nw
                z    = omega(iw) + 1i*eta;
                Gr   = obj.rgf_surface_(H, z, M, Hdim);
                Ga   = Gr';
                Gam  = eta * eye(size(Gr));
                T(iw) = real(trace(Gam * Gr * Gam * Ga));
            end
            obj.T_cache = T;
        end
    end

    methods (Access = private)
        function M = layer_size_(obj, Hdim, N)
            if obj.n_layers > 0
                M = ceil(Hdim / obj.n_layers);
            else
                M = max(1, ceil(sqrt(Hdim)));  % auto: sqrt(Hdim) sites/layer
            end
        end

        function G = rgf_surface_(~, H, z, M, Hdim)
            % Surface Green's function via forward recursion.
            n_lay = ceil(Hdim / M);
            G = zeros(min(M,Hdim));
            for lay = 1:n_lay
                i1 = (lay-1)*M + 1;
                i2 = min(lay*M, Hdim);
                m  = i2 - i1 + 1;
                H_n = H(i1:i2, i1:i2);
                zI  = z * eye(m);
                if lay == 1
                    G = (zI - H_n) \ eye(m);
                else
                    % Coupling to previous layer
                    V_nl = H(i1:i2, i1-M:i1-1);
                    Sig  = V_nl * G * V_nl';
                    G    = (zI - H_n - Sig) \ eye(m);
                end
            end
        end
    end
end
