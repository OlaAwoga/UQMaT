classdef Spectral
    % UQMaT.Observables.Conventional.Square2D.Spectral
    %   Spectral observables for the 2D square lattice.
    %
    %   Works with any UQMaT.Solver built from a
    %   UQMaT.Models.Lattice.Square2D.Hamiltonian system.
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params2D('Nx',20,'Ny',20,'t',1,'mu',-1);
    %     sys = UQMaT.System(p);
    %     slv = UQMaT.Solver.Diagonalization.Direct(sys);
    %     obs = UQMaT.Observables.Conventional.Square2D.Spectral(slv);
    %
    %     [w, D]         = obs.dos()
    %     [w, D_up,D_dn] = obs.dos_spin()         % requires spin=true
    %     [w, D_e, D_h]  = obs.dos_ph()           % requires nambu=true
    %     [w, L]         = obs.ldos()             % Nω × Ns  (all sites)
    %     [w, L]         = obs.ldos(sites)        % selected linear indices
    %     [w, Lmap]      = obs.ldos_map(omega_val) % Nx × Ny  at one energy
    %     [kd, Ek]       = obs.bands()            % band structure Γ→X→M→Γ
    %                       uses PBC k-space H(k) if bcx=bcy=pbc
    %     F              = obs.pair_amplitude()   % requires nambu=true
    %     F2D            = obs.pair_amplitude_map()  % Nx × Ny reshape
    %
    %   Notes
    %   -----
    %   - ldos_map: reshapes site index back to (ix,iy) grid for plotting.
    %   - bands: computed analytically from ε(k) = -2t(cos kx + cos ky) - μ
    %     when normal state; or by diagonalising H(k) for spinful/BdG.

    properties (SetAccess = private)
        solver
        p       % UQMaT.Params2D
        Nx
        Ny
        Ns
    end

    methods
        function obj = Spectral(solver)
            if ~isa(solver,'UQMaT.Solver.Base')
                error('UQMaT:Square2D:Spectral:badInput', ...
                    'Argument must be a UQMaT.Solver.Base subclass.');
            end
            if ~isa(solver.sys.params,'UQMaT.Params2D')
                error('UQMaT:Square2D:Spectral:notSquare2D', ...
                    'System must be built from UQMaT.Params2D.');
            end
            obj.solver = solver;
            obj.p      = solver.sys.params;
            obj.Nx     = obj.p.Nx;
            obj.Ny     = obj.p.Ny;
            obj.Ns     = obj.p.Ns;
        end

        % ================================================================
        % DOS (total, per physical site)
        % ================================================================
        function [omega, D] = dos(obj, varargin)
            [omega, eta] = obj.omegaGrid(varargin{:});
            [E, ~] = obj.solver.eigs();
            Nw = numel(omega);
            D  = zeros(Nw,1);
            for n = 1:numel(E)
                D = D + lorentz(omega, E(n), eta);
            end
            D = D / obj.Ns;
        end

        % ================================================================
        % Spin-resolved DOS (requires spin=true)
        % ================================================================
        function [omega, D_up, D_dn] = dos_spin(obj, varargin)
            obj.requireSpin();
            [omega, eta] = obj.omegaGrid(varargin{:});
            [E, V] = obj.solver.eigs();
            Ns = obj.Ns;
            Nw = numel(omega);
            D_up = zeros(Nw,1);
            D_dn = zeros(Nw,1);
            for n = 1:numel(E)
                lz   = lorentz(omega, E(n), eta);
                w_up = sum(abs(V(1:Ns,     n)).^2);
                w_dn = sum(abs(V(Ns+1:2*Ns, n)).^2);
                D_up = D_up + w_up * lz;
                D_dn = D_dn + w_dn * lz;
            end
            D_up = D_up / Ns;
            D_dn = D_dn / Ns;
        end

        % ================================================================
        % Particle-hole DOS (requires nambu=true)
        % ================================================================
        function [omega, D_e, D_h] = dos_ph(obj, varargin)
            obj.requireNambu();
            [omega, eta] = obj.omegaGrid(varargin{:});
            [E, V] = obj.solver.eigs();
            Ns = obj.Ns;
            Nw = numel(omega);
            D_e = zeros(Nw,1);
            D_h = zeros(Nw,1);
            for n = 1:numel(E)
                lz  = lorentz(omega, E(n), eta);
                w_e = sum(abs(V(1:2*Ns,      n)).^2);
                w_h = sum(abs(V(2*Ns+1:4*Ns, n)).^2);
                D_e = D_e + w_e * lz;
                D_h = D_h + w_h * lz;
            end
            D_e = D_e / Ns;
            D_h = D_h / Ns;
        end

        % ================================================================
        % LDOS  (Nω × Nsites selected)
        % ================================================================
        function [omega, L] = ldos(obj, sites, varargin)
            if nargin < 2 || isempty(sites), sites = 1:obj.Ns; end
            [omega, eta] = obj.omegaGrid(varargin{:});
            [E, V] = obj.solver.eigs();
            Ns = obj.Ns;

            % Electron block: rows 1:Ns (normal/spin-up) or 1:2Ns (spinful)
            if obj.p.nambu
                Vel = V(1:2*Ns, :);   % electron block (↑ + ↓)
                wt_fn = @(si,n) abs(Vel(si,n))^2 + abs(Vel(si+Ns,n))^2;
            elseif obj.p.spin
                Vel   = V(1:2*Ns, :);
                wt_fn = @(si,n) abs(Vel(si,n))^2 + abs(Vel(si+Ns,n))^2;
            else
                Vel   = V(1:Ns, :);
                wt_fn = @(si,n) abs(Vel(si,n))^2;
            end

            Nw = numel(omega);
            Nsel = numel(sites);
            L  = zeros(Nw, Nsel);
            for n = 1:numel(E)
                lz = lorentz(omega, E(n), eta);
                for si = 1:Nsel
                    i = sites(si);
                    L(:,si) = L(:,si) + wt_fn(i,n) * lz;
                end
            end
        end

        % ================================================================
        % LDOS map  at a single energy: returns Nx × Ny matrix
        % ================================================================
        function Lmap = ldos_map(obj, omega_val)
            [~, L] = obj.ldos(1:obj.Ns, 'omega', omega_val);
            Lmap   = reshape(L, obj.Ny, obj.Nx)';   % Nx × Ny
        end

        % ================================================================
        % Band structure along Γ → X → M → Γ
        % Extracted via Fourier transform of real-space H — no H(k) written.
        % Requires bcx='pbc', bcy='pbc'.
        % ================================================================
        function [kdist, Ek, klabels, kticks] = bands(obj, Nk)
            if nargin < 2, Nk = 200; end
            [kdist, Ek, klabels, kticks] = ...
                UQMaT.Models.Lattice.BoundaryConditions.BC.bands2D( ...
                    obj.solver.sys, Nk);
        end
            p = obj.p;

            [kpath, klabels, kticks] = ...
                UQMaT.Models.Lattice.Square2D.Geometry.kpath_square(Nk);
            Nkpts = size(kpath,1);
            t  = p.t;  t2 = p.t2;  mu = p.mu;

            % Band count
            if p.nambu
                Nbands = 4;
            elseif p.spin
                Nbands = 2;
            else
                Nbands = 1;
            end

            Ek    = zeros(Nkpts, Nbands);
            kdist = zeros(Nkpts,1);
            for ik = 2:Nkpts
                kdist(ik) = kdist(ik-1) + norm(kpath(ik,:) - kpath(ik-1,:));
            end

            for ik = 1:Nkpts
                kx = kpath(ik,1);  ky = kpath(ik,2);
                ek = -2*t*(cos(kx)+cos(ky)) - 4*t2*cos(kx)*cos(ky) - mu;

                if ~p.spin && ~p.nambu
                    Ek(ik,1) = ek;

                elseif p.spin && ~p.nambu
                    sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
                    HkN = ek*eye(2) + p.Bz*sig.sz + p.Bx*sig.sx + p.By*sig.sy;
                    if p.lambda_R ~= 0
                        % Rashba: λ_R(sin(ky) σ_x - sin(kx) σ_y)
                        HkN = HkN + p.lambda_R * ...
                            (sin(ky)*sig.sx - sin(kx)*sig.sy);
                    end
                    Ek(ik,:) = sort(real(eig(HkN)))';

                else
                    % Nambu BdG H(k): 4×4
                    sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
                    sy  = sig.sy;  sx = sig.sx;  sz = sig.sz;
                    HkN = ek*eye(2) + p.Bz*sz + p.Bx*sx + p.By*sy;
                    if p.lambda_R ~= 0
                        HkN = HkN + p.lambda_R*(sin(ky)*sx - sin(kx)*sy);
                    end
                    % Singlet pairing at k
                    Ds = mean(p.Delta_s_vec);
                    Dk = Ds * (1i*sy);  % s-wave: k-independent
                    nb = UQMaT.Models.Lattice.Nambu.Basis;
                    Hk = nb.assemble(HkN, Dk);
                    Ek(ik,:) = sort(real(eig(Hk)))';
                end
            end
        end

        % ================================================================
        % Pair amplitude  F(i) = <c_{i↓} c_{i↑}>  (requires nambu=true)
        % ================================================================
        function F = pair_amplitude(obj)
            obj.requireNambu();
            [E, V] = obj.solver.eigs();
            Ns  = obj.Ns;
            F   = zeros(Ns, 1);
            occ = find(E < 0);
            for n = occ'
                u_up = V(1:Ns,       n);   % electron ↑
                v_dn = V(3*Ns+1:4*Ns, n);  % hole ↓
                F = F + conj(u_up) .* v_dn;
            end
        end

        % ================================================================
        % Pair amplitude as Nx × Ny map
        % ================================================================
        function F2D = pair_amplitude_map(obj)
            F   = obj.pair_amplitude();
            F2D = reshape(F, obj.Ny, obj.Nx)';   % Nx × Ny
        end

        % ================================================================
        % Spectrum
        % ================================================================
        function [E, V] = spectrum(obj)
            [E, V] = obj.solver.eigs();
        end

        % ================================================================
        % Display
        % ================================================================
        function disp(obj)
            p = obj.p;
            fprintf('  UQMaT.Observables.Conventional.Square2D.Spectral\n');
            fprintf('  Grid: %d×%d  pairing=%s  spin=%d  nambu=%d\n', ...
                p.Nx, p.Ny, p.pairing, p.spin, p.nambu);
        end
    end

    % ---- Private helpers -----------------------------------------------
    methods (Access = private)
        function [omega, eta] = omegaGrid(obj, varargin)
            p   = obj.p;
            eta = p.eta;
            idx = find(cellfun(@(x) ischar(x) && strcmpi(x,'omega'), varargin));
            if ~isempty(idx)
                omega = varargin{idx+1}(:);
                return
            end
            E   = obj.solver.sys.E;
            pad = max(0.2*(E(end)-E(1)), 4*eta);
            omega = linspace(E(1)-pad, E(end)+pad, p.Nomega)';
        end

        function requireSpin(obj)
            if ~obj.p.spin
                error('UQMaT:Square2D:Spectral:noSpin', ...
                    'This method requires spin=true in Params2D.');
            end
        end

        function requireNambu(obj)
            if ~obj.p.nambu
                error('UQMaT:Square2D:Spectral:noNambu', ...
                    'This method requires nambu=true (Nambu BdG) in Params2D.');
            end
        end
    end
end

function y = lorentz(omega, E0, eta)
    y = (eta/pi) ./ ((omega - E0).^2 + eta^2);
end
