classdef Eilenberger
    % UQMaT.Models.Quasiclassical.Eilenberger
    %   Eilenberger equation solver for clean-limit quasiclassical SC.
    %
    %   The Eilenberger equation describes quasiclassical Green's functions
    %   in the clean limit (mean free path >> coherence length):
    %
    %     ℏv_F · ∇ĝ + [Ĥ_qc, ĝ] = 0
    %
    %   where ĝ is the 2×2 quasiclassical Green's function matrix:
    %     ĝ = [ g,   f  ]
    %         [-f̃, -g̃ ]
    %   with g = normal GF, f = anomalous GF (pair amplitude)
    %
    %   Normalisation constraint: ĝ² = -π²·Î
    %
    %   This implementation uses the Ricatti parametrisation:
    %     ĝ = -iπN_F [ (1-γγ̃)⁻¹(1+γγ̃),  2(1-γγ̃)⁻¹γ  ]
    %                 [ -2(1-γ̃γ)⁻¹γ̃,     -(1-γ̃γ)⁻¹(1+γ̃γ) ]
    %
    %   Reference: Eilenberger G., Z. Phys. 214, 195 (1968)
    %              Schopohl N., Maki K., PRB 52, 490 (1995) (Ricatti)
    %
    %   Usage:
    %     el = UQMaT.Models.Quasiclassical.Eilenberger(params)
    %     el = el.solve()           % solve Eilenberger equation
    %     [g, f] = el.green(x)     % Green's function at position x
    %     Delta  = el.gap_profile() % self-consistent gap
    %
    %   Params fields:
    %     qc_vF      : Fermi velocity (default 1.0)
    %     qc_Delta0  : bulk gap (sets initial profile)
    %     qc_lambda  : pairing interaction strength (for self-consistency)
    %     qc_geometry: 'bulk' | 'slab' | 'vortex'
    %     qc_Lx      : system length in x (for slab)
    %     qc_Nx      : number of discretisation points
    %     qc_T       : temperature (default 0 = zero temperature)
    %     qc_omega_n : Matsubara frequencies (auto-generated if empty)
    %     qc_N_omega : number of Matsubara frequencies (default 50)

    %  Citation
    %  --------
    %  Eilenberger G., Z. Phys. 214, 195 (1968). [Eilenberger equation]
    %  Serene J.W., Rainer D., Phys. Rep. 101, 221 (1983). [quasiclassical review]
    %  Schopohl N., Maki K., PRB 52, 490 (1995). [Ricatti parametrisation]
    %  Belzig W. et al., Superlattices Microstruct. 25, 1251 (1999). [review]
    properties
        params          % UQMaT.Params
        x      (:,1) double   % spatial grid
        Delta  (:,1) double   % gap profile
        gamma  (:,:) double   % Ricatti amplitude (Nx × N_omega)
        gamma_tilde (:,:) double
        converged (1,1) logical = false
    end

    methods

        function obj = Eilenberger(p)
            obj.params = p;
            obj = obj.init_grid();
            obj = obj.init_delta();
        end

        function obj = solve(obj, varargin)
            % Solve Eilenberger equation via Ricatti parametrisation.
            %   obj = el.solve()
            %   obj = el.solve('tol',1e-6,'maxiter',100,'verbose',true)
            opt = struct('tol',1e-6,'maxiter',100,'verbose',false,'selfconsistent',false);
            for k=1:2:numel(varargin), opt.(varargin{k})=varargin{k+1}; end
            p = obj.params;

            omega_n = obj.matsubara_freqs();
            Nx      = numel(obj.x);
            N_om    = numel(omega_n);

            obj.gamma       = zeros(Nx, N_om);
            obj.gamma_tilde = zeros(Nx, N_om);

            for sc_iter = 1:opt.maxiter
                Delta_old = obj.Delta;

                % Solve Ricatti ODE for each Matsubara frequency
                for n = 1:N_om
                    [obj.gamma(:,n), obj.gamma_tilde(:,n)] = ...
                        obj.ricatti_solve(omega_n(n), obj.Delta, p);
                end

                % Update gap if self-consistent
                if opt.selfconsistent
                    obj = obj.update_gap(omega_n);
                    err = max(abs(obj.Delta - Delta_old));
                    if opt.verbose
                        fprintf('  Eilenberger iter %d: max|ΔΔ|=%.3e\n', sc_iter, err);
                    end
                    if err < opt.tol
                        obj.converged = true; break;
                    end
                else
                    obj.converged = true; break;
                end
            end
        end

        function [g, f] = green(obj, x_query)
            % Return quasiclassical Green's function at position x.
            % g: N_omega×1 normal GF;  f: N_omega×1 anomalous GF
            omega_n = obj.matsubara_freqs();
            ix = UQMaT.Models.Quasiclassical.Eilenberger.nearest_idx(x_query, obj.x);
            gam   = obj.gamma(ix,:)';
            gam_t = obj.gamma_tilde(ix,:)';
            denom = 1 + gam .* gam_t;
            g = -1i*pi * (1 - gam.*gam_t) ./ denom;
            f = -1i*pi * 2 * gam ./ denom;
        end

        function Delta = gap_profile(obj)
            Delta = obj.Delta;
        end

        function [g, f, x] = profile(obj)
            % Full Green's function profile.
            x = obj.x;
            omega_n = obj.matsubara_freqs();
            N_om = numel(omega_n);
            Nx   = numel(x);
            g    = zeros(Nx, N_om);
            f    = zeros(Nx, N_om);
            for ix = 1:Nx
                [g(ix,:), f(ix,:)] = obj.green(x(ix));
            end
        end

        function disp(obj)
            p = obj.params;
            fprintf('\n  UQMaT Eilenberger solver\n');
            fprintf('  geometry : %s\n', obj.get_field(p,'qc_geometry','slab'));
            fprintf('  Nx       : %d\n', numel(obj.x));
            fprintf('  vF       : %.4g\n', obj.get_field(p,'qc_vF',1));
            fprintf('  Delta0   : %.4g\n', obj.get_field(p,'qc_Delta0',0.1));
            fprintf('  T        : %.4g\n', obj.get_field(p,'qc_T',0));
            fprintf('  converged: %d\n\n', obj.converged);
        end

    end

    methods (Access = private)

        function obj = init_grid(obj)
            p  = obj.params;
            Nx = obj.get_field(p,'qc_Nx',100);
            Lx = obj.get_field(p,'qc_Lx',10);
            obj.x = linspace(0, Lx, Nx)';
        end

        function obj = init_delta(obj)
            p      = obj.params;
            Delta0 = obj.get_field(p,'qc_Delta0',0.1);
            geom   = obj.get_field(p,'qc_geometry','bulk');
            Nx     = numel(obj.x);

            switch geom
                case 'bulk'
                    obj.Delta = Delta0 * ones(Nx,1);
                case 'slab'
                    % Step profile: full gap everywhere initially
                    obj.Delta = Delta0 * ones(Nx,1);
                case 'vortex'
                    % Radial profile: Delta(r) = Delta0 * tanh(r/xi)
                    r = obj.x;
                    xi = 1;   % coherence length in units of vF/Delta0
                    obj.Delta = Delta0 * tanh(r/xi);
                otherwise
                    obj.Delta = Delta0 * ones(Nx,1);
            end
        end

        function omega_n = matsubara_freqs(obj)
            p     = obj.params;
            T     = obj.get_field(p,'qc_T',0.01);   % avoid T=0 singularity
            N_om  = obj.get_field(p,'qc_N_omega',50);
            if isfield(p,'qc_omega_n') && ~isempty(p.qc_omega_n)
                omega_n = p.qc_omega_n(:);
            else
                n       = (0:N_om-1)';
                omega_n = pi*T*(2*n+1);  % fermionic Matsubara: (2n+1)πT
            end
        end

        function [gam, gam_t] = ricatti_solve(obj, omega_n, Delta, p)
            % Solve Ricatti ODE: ℏvF * d(γ)/dx = 2ω_n*γ + Δ̃ - Δ*γ²
            % on spatial grid x using RK4.
            vF  = obj.get_field(p,'qc_vF',1);
            Nx  = numel(obj.x);
            dx  = obj.x(2) - obj.x(1);

            gam   = zeros(Nx,1);
            gam_t = zeros(Nx,1);

            % Initial conditions (bulk BCS)
            Delta0_bulk = obj.get_field(p,'qc_Delta0',0.1);
            denom_bulk  = sqrt(omega_n^2 + Delta0_bulk^2);
            gam(1)      = Delta0_bulk / (omega_n + denom_bulk);
            gam_t(Nx)   = conj(gam(1));

            % Forward integration for gamma
            for ix = 1:Nx-1
                g = gam(ix); D = Delta(ix); Dt = conj(Delta(ix));
                rhs = (2*omega_n*g + Dt - D*g^2) / vF;
                gam(ix+1) = gam(ix) + dx*rhs;
            end
            % Backward integration for gamma_tilde
            for ix = Nx:-1:2
                gt = gam_t(ix); D = Delta(ix); Dt = conj(Delta(ix));
                rhs = -(2*omega_n*gt + D - Dt*gt^2) / vF;
                gam_t(ix-1) = gam_t(ix) - dx*rhs;
            end
        end

        function obj = update_gap(obj, omega_n)
            p      = obj.params;
            lambda = obj.get_field(p,'qc_lambda',1.0);
            T      = obj.get_field(p,'qc_T',0.01);
            Nx     = numel(obj.x);
            N_om   = numel(omega_n);

            Delta_new = zeros(Nx,1);
            for ix = 1:Nx
                s = 0;
                for n = 1:N_om
                    gam   = obj.gamma(ix,n);
                    gam_t = obj.gamma_tilde(ix,n);
                    f_n   = -1i*pi*2*gam/(1 + gam*gam_t);
                    s     = s + real(f_n);
                end
                Delta_new(ix) = lambda * T * s;
            end
            obj.Delta = 0.5*obj.Delta + 0.5*Delta_new;  % linear mixing
        end

    end

    methods (Static, Access=private)
        function v = get_field(p, fname, default)
            if isfield(p, fname), v = p.(fname);
            else, v = default;
            end
        end
        function ix = nearest_idx(x_query, x_grid)
            [~,ix] = min(abs(x_grid - x_query));
        end
    end
end
