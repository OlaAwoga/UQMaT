classdef Nambu
    % UQMaT.Observables.Conventional.Nambu
    %   Spectral observables for superconducting systems in the Nambu basis.
    %
    %   The Nambu basis Ψ = (c↑, c↓, c↑†, c↓†)ᵀ IS the superconducting
    %   basis — this class is the standard observable class for any system
    %   built with sc=true.
    %
    %   Row blocks (4N × 4N):
    %     1    : N     →  electron ↑
    %     N+1  : 2N    →  electron ↓
    %     2N+1 : 3N    →  hole ↑
    %     3N+1 : 4N    →  hole ↓
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params('N',80,'sc',true,'Delta_s',0.3,'Bz',0.2);
    %     sys = UQMaT.System(p);
    %     slv = UQMaT.Solver.Diagonalization.Direct(sys);
    %     obs = UQMaT.Observables.Conventional.Nambu(slv);
    %
    %     [w, D]           = obs.dos()
    %     [w, D_e, D_h]    = obs.dos_ph()        % particle / hole DOS
    %     [w, D_up, D_dn]  = obs.dos_spin()      % spin ↑ / ↓ DOS
    %     [w, L]           = obs.ldos()          % total LDOS, Nω × N
    %     [w, L]           = obs.ldos(sites)
    %     [w, L_up, L_dn]  = obs.ldos_spin(sites)
    %     [w, L_e, L_h]    = obs.ldos_ph(sites)
    %     [w, Sx, Sy, Sz]  = obs.spin_texture(sites)
    %     F                = obs.pair_amplitude() % <c↑c↓> per site, N×1
    %     [E, V]           = obs.spectrum()

    properties (SetAccess = private)
        solver
        N      (1,1) double
    end

    methods
        function obj = Nambu(solver)
            if ~isa(solver,'UQMaT.Solver.Base')
                error('UQMaT:Observables:Nambu:badInput', ...
                    'Argument must be a UQMaT.Solver.Base subclass.');
            end
            if ~solver.sys.params.bdg
                error('UQMaT:Observables:Nambu:notSC', ...
                    'System must be built with sc=true.');
            end
            obj.solver = solver;
            obj.N      = solver.sys.params.N;
        end

        % ================================================================
        % DOS — total, per physical site
        % ================================================================
        function [omega, D] = dos(obj, varargin)
            [omega, eta] = obj.grid(varargin{:});
            E  = obj.solver.eigvals();
            D  = zeros(numel(omega),1);
            for n = 1:numel(E)
                D = D + lor(omega, E(n), eta);
            end
            D = D / obj.N;
        end

        % ================================================================
        % Particle / hole resolved DOS
        % ================================================================
        function [omega, D_e, D_h] = dos_ph(obj, varargin)
            [omega, eta] = obj.grid(varargin{:});
            [E, V] = obj.solver.eigs();
            N  = obj.N;
            nb = UQMaT.Models.Lattice.Nambu.Basis;
            Nw = numel(omega);
            D_e = zeros(Nw,1);  D_h = zeros(Nw,1);
            for n = 1:numel(E)
                lz  = lor(omega, E(n), eta);
                w_e = sum(abs(V(nb.elec_block(N), n)).^2);
                w_h = sum(abs(V(nb.hole_block(N), n)).^2);
                D_e = D_e + w_e*lz;
                D_h = D_h + w_h*lz;
            end
            D_e = D_e/N;  D_h = D_h/N;
        end

        % ================================================================
        % Spin-resolved DOS
        % ================================================================
        function [omega, D_up, D_dn] = dos_spin(obj, varargin)
            [omega, eta] = obj.grid(varargin{:});
            [E, V] = obj.solver.eigs();
            N  = obj.N;
            nb = UQMaT.Models.Lattice.Nambu.Basis;
            Nw = numel(omega);
            D_up = zeros(Nw,1);  D_dn = zeros(Nw,1);
            for n = 1:numel(E)
                lz     = lor(omega, E(n), eta);
                w_up   = sum(abs(V(nb.elec_up(N), n)).^2);
                w_dn   = sum(abs(V(nb.elec_dn(N), n)).^2);
                D_up   = D_up + w_up*lz;
                D_dn   = D_dn + w_dn*lz;
            end
            D_up = D_up/N;  D_dn = D_dn/N;
        end

        % ================================================================
        % LDOS — electron spectral weight per site
        % ================================================================
        function [omega, L] = ldos(obj, sites, varargin)
            if nargin < 2 || isempty(sites), sites = 1:obj.N; end
            [omega, eta] = obj.grid(varargin{:});
            [E, V] = obj.solver.eigs();
            N  = obj.N;
            Nw = numel(omega);
            Ns = numel(sites);
            L  = zeros(Nw, Ns);
            for n = 1:numel(E)
                lz = lor(omega, E(n), eta);
                for si = 1:Ns
                    i  = sites(si);
                    % Sum electron ↑ and ↓ weight at site i
                    w  = abs(V(i, n))^2 + abs(V(i+N, n))^2;
                    L(:,si) = L(:,si) + w*lz;
                end
            end
        end

        % ================================================================
        % Spin-resolved LDOS
        % ================================================================
        function [omega, L_up, L_dn] = ldos_spin(obj, sites, varargin)
            if nargin < 2 || isempty(sites), sites = 1:obj.N; end
            [omega, eta] = obj.grid(varargin{:});
            [E, V] = obj.solver.eigs();
            N  = obj.N;
            Nw = numel(omega);
            Ns = numel(sites);
            L_up = zeros(Nw,Ns);  L_dn = zeros(Nw,Ns);
            for n = 1:numel(E)
                lz = lor(omega, E(n), eta);
                for si = 1:Ns
                    i = sites(si);
                    L_up(:,si) = L_up(:,si) + abs(V(i,   n))^2 * lz;
                    L_dn(:,si) = L_dn(:,si) + abs(V(i+N, n))^2 * lz;
                end
            end
        end

        % ================================================================
        % Particle-hole LDOS
        % ================================================================
        function [omega, L_e, L_h] = ldos_ph(obj, sites, varargin)
            if nargin < 2 || isempty(sites), sites = 1:obj.N; end
            [omega, eta] = obj.grid(varargin{:});
            [E, V] = obj.solver.eigs();
            N  = obj.N;
            Nw = numel(omega);
            Ns = numel(sites);
            L_e = zeros(Nw,Ns);  L_h = zeros(Nw,Ns);
            for n = 1:numel(E)
                lz = lor(omega, E(n), eta);
                for si = 1:Ns
                    i  = sites(si);
                    we = abs(V(i,     n))^2 + abs(V(i+N,   n))^2;   % e↑+e↓
                    wh = abs(V(i+2*N, n))^2 + abs(V(i+3*N, n))^2;   % h↑+h↓
                    L_e(:,si) = L_e(:,si) + we*lz;
                    L_h(:,si) = L_h(:,si) + wh*lz;
                end
            end
        end

        % ================================================================
        % Spin texture  <S_α>(ω, site)
        % ================================================================
        function [omega, Sx, Sy, Sz] = spin_texture(obj, sites, varargin)
            if nargin < 2 || isempty(sites), sites = 1:obj.N; end
            [omega, eta] = obj.grid(varargin{:});
            [E, V] = obj.solver.eigs();
            N   = obj.N;
            sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
            Nw  = numel(omega);
            Ns  = numel(sites);
            Sx  = zeros(Nw,Ns);
            Sy  = zeros(Nw,Ns);
            Sz  = zeros(Nw,Ns);
            for n = 1:numel(E)
                lz = lor(omega, E(n), eta);
                for si = 1:Ns
                    i   = sites(si);
                    psi = [V(i,n); V(i+N,n)];   % [↑, ↓] at site i
                    Sx(:,si) = Sx(:,si) + real(psi'*sig.sx*psi)*lz;
                    Sy(:,si) = Sy(:,si) + real(psi'*sig.sy*psi)*lz;
                    Sz(:,si) = Sz(:,si) + real(psi'*sig.sz*psi)*lz;
                end
            end
        end

        % ================================================================
        % Pair amplitude  F(i) = <c_{i↓} c_{i↑}>
        %
        % From BdG eigenstates:
        %   F(i) = Σ_{E_n < 0}  u_{n↑}(i) · v_{n↓}*(i)
        %
        % where u_{n↑} = V(i, n)     (electron ↑ component)
        %       v_{n↓} = V(i+3N, n)  (hole ↓ component, rows 3N+1:4N)
        %
        % Bug fix: previous code had conj(u)*v; correct is u*conj(v).
        % ================================================================
        function F = pair_amplitude(obj)
            [E, V] = obj.solver.eigs();
            N   = obj.N;
            F   = zeros(N,1);
            occ = find(E < 0);
            for n = occ'
                u_up = V(1:N,       n);   % electron ↑,  rows 1:N
                v_dn = V(3*N+1:4*N, n);   % hole ↓,      rows 3N+1:4N
                % F(i) = u↑(i) · v↓*(i)   [fixed sign]
                F = F + u_up .* conj(v_dn);
            end
        end

        % ================================================================
        % Spectrum
        % ================================================================
        function [E, V] = spectrum(obj)
            [E, V] = obj.solver.eigs();
        end

        function disp(obj)
            p = obj.solver.sys.params;
            fprintf('  UQMaT.Observables.Conventional.Nambu\n');
            fprintf('  Basis : Ψ = (c↑, c↓, c↑†, c↓†)ᵀ  (sc=true)\n');
            fprintf('  N=%d, pairing=%s\n', p.N, p.pairing_mode);
        end
    end

    methods (Access = private)
        function [omega, eta] = grid(obj, varargin)
            p   = obj.solver.sys.params;
            eta = p.eta;
            idx = find(cellfun(@(x) ischar(x) && strcmpi(x,'omega'), varargin));
            if ~isempty(idx)
                omega = varargin{idx+1}(:); return
            end
            E   = obj.solver.sys.E;
            pad = max(0.2*(E(end)-E(1)), 4*eta);
            omega = linspace(E(1)-pad, E(end)+pad, 500)';
        end
    end
end

function y = lor(omega, E0, eta)
    y = (eta/pi) ./ ((omega-E0).^2 + eta^2);
end
