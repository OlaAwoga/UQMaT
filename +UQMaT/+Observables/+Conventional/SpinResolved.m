classdef SpinResolved
    % UQMaT.Observables.Conventional.SpinResolved
    %   Spin-resolved spectral quantities for spinful systems (v1.4).
    %   Requires a system built with spin=true (Zeeman or Rashba active).
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params('spin',true,'Bz',0.5,'lambda_R',0.3,'N',60);
    %     sys = UQMaT.System(p);
    %     slv = UQMaT.Solver.Diagonalization.Direct(sys);
    %     obs = UQMaT.Observables.Conventional.SpinResolved(slv);
    %
    %     [w, D_up, D_dn]   = obs.dos_spin()     % spin-↑ and ↓ DOS
    %     [w, L_up, L_dn]   = obs.ldos_spin()    % spin-resolved LDOS
    %     [w, L_up, L_dn]   = obs.ldos_spin(sites)
    %     [w, Sx, Sy, Sz]   = obs.spin_texture()  % <S_alpha>(omega) per site
    %
    %   Notes
    %   -----
    %   - Spin-↑ block: rows 1:N,      Spin-↓ block: rows N+1:2N.
    %   - For BdG (sc=true): electron block is 1:2N (↑+↓ particle), hole 2N+1:4N.
    %   - spin_texture returns Nomega × Nsites for each component.

    properties (SetAccess = private)
        solver
    end

    methods
        function obj = SpinResolved(solver)
            if ~isa(solver,'UQMaT.Solver.Base')
                error('UQMaT:SpinResolved:badInput', ...
                    'Argument must be a UQMaT.Solver.Base subclass.');
            end
            if ~solver.sys.params.spin
                error('UQMaT:SpinResolved:noSpin', ...
                    'System must be built with spin=true.');
            end
            obj.solver = solver;
        end

        % ---- Spin-resolved DOS -----------------------------------------
        function [omega, D_up, D_dn] = dos_spin(obj, varargin)
            [omega, eta] = obj.omegaGrid(varargin{:});
            [E, V] = obj.solver.eigs();
            p   = obj.solver.sys.params;
            N   = p.N;
            Nw  = numel(omega);

            D_up = zeros(Nw,1);
            D_dn = zeros(Nw,1);

            for n = 1:numel(E)
                lz = lorentz(omega, E(n), eta);
                wt_up = sum(abs(V(1:N, n)).^2);
                wt_dn = sum(abs(V(N+1:2*N, n)).^2);
                D_up = D_up + wt_up * lz;
                D_dn = D_dn + wt_dn * lz;
            end

            D_up = D_up / N;
            D_dn = D_dn / N;
        end

        % ---- Spin-resolved LDOS ----------------------------------------
        function [omega, L_up, L_dn] = ldos_spin(obj, sites, varargin)
            if nargin < 2 || isempty(sites)
                sites = 1:obj.solver.sys.params.N;
            end
            [omega, eta] = obj.omegaGrid(varargin{:});
            [E, V] = obj.solver.eigs();
            p  = obj.solver.sys.params;
            N  = p.N;
            Nw = numel(omega);
            Ns = numel(sites);

            L_up = zeros(Nw, Ns);
            L_dn = zeros(Nw, Ns);

            for n = 1:numel(E)
                lz = lorentz(omega, E(n), eta);
                for si = 1:Ns
                    i = sites(si);
                    w_up = abs(V(i,   n))^2;
                    w_dn = abs(V(i+N, n))^2;
                    L_up(:,si) = L_up(:,si) + w_up * lz;
                    L_dn(:,si) = L_dn(:,si) + w_dn * lz;
                end
            end
        end

        % ---- Spin texture <S_alpha>(omega, site) -----------------------
        function [omega, Sx, Sy, Sz] = spin_texture(obj, sites, varargin)
            % Spin expectation value per energy and site:
            %   S_alpha(omega, i) = sum_n |<n|S_alpha|i>|^2 * delta(omega-E_n)
            if nargin < 2 || isempty(sites)
                sites = 1:obj.solver.sys.params.N;
            end
            [omega, eta] = obj.omegaGrid(varargin{:});
            [E, V] = obj.solver.eigs();
            N  = obj.solver.sys.params.N;
            Nw = numel(omega);
            Ns = numel(sites);

            Sx = zeros(Nw, Ns);
            Sy = zeros(Nw, Ns);
            Sz = zeros(Nw, Ns);

            % Pauli matrices (site-local): sigma acts on [↑_i, ↓_i] subspace
            for n = 1:numel(E)
                lz = lorentz(omega, E(n), eta);
                for si = 1:Ns
                    i = sites(si);
                    psi_up = V(i,   n);
                    psi_dn = V(i+N, n);
                    spinor = [psi_up; psi_dn];

                    sx = real(spinor' * [0 1; 1 0]  * spinor);
                    sy = real(spinor' * [0 -1i; 1i 0] * spinor);
                    sz = real(spinor' * [1 0; 0 -1] * spinor);

                    Sx(:,si) = Sx(:,si) + sx * lz;
                    Sy(:,si) = Sy(:,si) + sy * lz;
                    Sz(:,si) = Sz(:,si) + sz * lz;
                end
            end
        end

        function disp(obj)
            p = obj.solver.sys.params;
            fprintf('  UQMaT.Observables.Conventional.SpinResolved\n');
            fprintf('  Bx=%.4g  By=%.4g  Bz=%.4g  lambda_R=%.4g\n', ...
                p.B(1), p.B(2), p.B(3), p.lambda_R);
        end
    end

    methods (Access = private)
        function [omega, eta] = omegaGrid(obj, varargin)
            p   = obj.solver.sys.params;
            eta = p.eta;
            idx = find(cellfun(@(x) ischar(x) && strcmpi(x,'omega'), varargin));
            if ~isempty(idx)
                omega = varargin{idx+1}(:);
                return
            end
            E   = obj.solver.sys.E;
            pad = max(0.2*(E(end)-E(1)), 4*eta);
            omega = linspace(E(1)-pad, E(end)+pad, 500)';
        end
    end
end

function y = lorentz(omega, E0, eta)
    y = (eta/pi) ./ ((omega - E0).^2 + eta^2);
end
