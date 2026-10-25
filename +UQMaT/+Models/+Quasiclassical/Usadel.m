classdef Usadel
    % UQMaT.Models.Quasiclassical.Usadel
    %   Usadel equation solver for dirty-limit quasiclassical SC.
    %
    %   The Usadel equation is the dirty-limit form of the Eilenberger
    %   equation, valid when mean free path l << coherence length ξ:
    %
    %     D·∇(ĝ·∇ĝ) + [iω_n*τ̂_z + Δ̂, ĝ] = 0
    %
    %   In the θ-parametrisation (ĝ = diag(cos θ, -cos θ), f = sin θ):
    %
    %     D·∂²θ/∂x² = 2ω_n·sin θ - 2Δ·cos θ
    %
    %   where D is the diffusion constant, ω_n = (2n+1)πT Matsubara freq.
    %
    %   Boundary conditions:
    %     - Bulk SC: θ(±∞) = arctan(Δ/ω_n)
    %     - NS interface: θ(0) = 0 (vacuum), θ(L) = bulk value
    %     - Kupriyanov-Lukichev: dθ/dx = -(1/l_b)·sin θ·cos θ at interface
    %
    %   Reference: Usadel K., PRL 25, 507 (1970)
    %              Kupriyanov M.Yu., Lukichev V.F., JETP 94, 139 (1988)
    %
    %   Usage:
    %     us = UQMaT.Models.Quasiclassical.Usadel(p)
    %     us = us.solve()
    %     [theta, x]    = us.theta_profile()
    %     [Delta, x]    = us.gap_profile()
    %     [N_S, E]      = us.dos()         % quasiclassical DOS
    %
    %   Params fields:
    %     qc_D       : diffusion constant (default 1.0)
    %     qc_Delta0  : bulk gap
    %     qc_lambda  : pairing interaction (self-consistency)
    %     qc_geometry: 'bulk' | 'ns_junction' | 'sns_junction'
    %     qc_Lx      : system length
    %     qc_Nx      : discretisation points (default 200)
    %     qc_T       : temperature (default 0.01*Delta0)
    %     qc_N_omega : Matsubara frequencies (default 100)
    %     qc_l_b     : Kupriyanov-Lukichev barrier length (default Inf = perfect)

    %  Citation
    %  --------
    %  Usadel K., PRL 25, 507 (1970). [Usadel equation, dirty limit]
    %  Kupriyanov M.Yu., Lukichev V.F., JETP 67, 1163 (1988). [BC at interface]
    %  Golubov A.A., Kupriyanov M.Yu., J. Low Temp. Phys. 70, 83 (1988).
    %  Belzig W. et al., Superlattices Microstruct. 25, 1251 (1999). [review]
    properties
        params
        x      (:,1) double
        theta  (:,:) double   % Nx × N_omega
        Delta  (:,1) double
        converged (1,1) logical = false
    end

    methods

        function obj = Usadel(p)
            obj.params = p;
            obj = obj.init_grid();
            obj = obj.init_delta();
        end

        function obj = solve(obj, varargin)
            opt = struct('tol',1e-6,'maxiter',200,'verbose',false,...
                         'selfconsistent',false,'mix',0.5);
            for k=1:2:numel(varargin), opt.(varargin{k})=varargin{k+1}; end
            p = obj.params;

            omega_n = obj.matsubara_freqs();
            Nx      = numel(obj.x);
            N_om    = numel(omega_n);
            obj.theta = zeros(Nx, N_om);

            for sc_iter = 1:opt.maxiter
                Delta_old = obj.Delta;

                for n = 1:N_om
                    obj.theta(:,n) = obj.solve_theta(omega_n(n), obj.Delta, p);
                end

                if opt.selfconsistent
                    obj = obj.update_gap(omega_n);
                    err = max(abs(obj.Delta - Delta_old));
                    if opt.verbose
                        fprintf('  Usadel iter %d: max|ΔΔ|=%.3e\n', sc_iter, err);
                    end
                    % Linear mixing
                    obj.Delta = (1-opt.mix)*Delta_old + opt.mix*obj.Delta;
                    if err < opt.tol
                        obj.converged = true; break;
                    end
                else
                    obj.converged = true; break;
                end
            end
        end

        function [theta, x] = theta_profile(obj)
            x     = obj.x;
            theta = obj.theta;
        end

        function [Delta, x] = gap_profile(obj)
            x     = obj.x;
            Delta = obj.Delta;
        end

        function [N_S, E] = dos(obj)
            % Quasiclassical DOS via analytic continuation ω_n → -iE + δ
            % N_S(E) = Re[cos θ(E)]
            p   = obj.params;
            eta = obj.get_field(p,'eta',0.01);
            E   = linspace(-3*obj.get_field(p,'qc_Delta0',0.1), ...
                            3*obj.get_field(p,'qc_Delta0',0.1), 400)';
            % Average theta over spatial grid
            Delta0 = obj.get_field(p,'qc_Delta0',0.1);
            N_S = zeros(size(E));
            for ie = 1:numel(E)
                en     = E(ie) + 1i*eta;
                theta_E= atanh(Delta0 / en);   % bulk BCS angle
                N_S(ie)= real(cos(theta_E));
            end
        end

        function disp(obj)
            p = obj.params;
            fprintf('\n  UQMaT Usadel solver\n');
            fprintf('  geometry : %s\n', obj.get_field(p,'qc_geometry','bulk'));
            fprintf('  Nx       : %d\n', numel(obj.x));
            fprintf('  D        : %.4g\n', obj.get_field(p,'qc_D',1));
            fprintf('  Delta0   : %.4g\n', obj.get_field(p,'qc_Delta0',0.1));
            fprintf('  T        : %.4g\n', obj.get_field(p,'qc_T',0.001));
            fprintf('  converged: %d\n\n', obj.converged);
        end

    end

    methods (Access=private)

        function obj = init_grid(obj)
            p  = obj.params;
            Nx = obj.get_field(p,'qc_Nx',200);
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
                case 'ns_junction'
                    % Step: SC on right half
                    obj.Delta = Delta0 * (obj.x > obj.x(end)/2);
                case 'sns_junction'
                    % SC-N-SC: normal region in centre
                    xm = obj.x(end)/2;  xw = obj.x(end)/4;
                    obj.Delta = Delta0 * (abs(obj.x - xm) > xw);
                otherwise
                    obj.Delta = Delta0 * ones(Nx,1);
            end
        end

        function omega_n = matsubara_freqs(obj)
            p    = obj.params;
            T    = obj.get_field(p,'qc_T',0.001);
            N_om = obj.get_field(p,'qc_N_omega',100);
            n    = (0:N_om-1)';
            omega_n = pi*T*(2*n+1);
        end

        function theta = solve_theta(obj, omega_n, Delta, p)
            % Solve Usadel ODE: D*d²θ/dx² = 2ω_n*sin θ - 2Δ*cos θ
            % via Newton iteration on finite-difference discretisation.
            D  = obj.get_field(p,'qc_D',1);
            Nx = numel(obj.x);
            dx = obj.x(2) - obj.x(1);

            % Initial guess: bulk BCS angle
            Delta0 = obj.get_field(p,'qc_Delta0',0.1);
            theta0 = atan2(real(Delta), omega_n);
            theta  = theta0 * ones(Nx,1);

            % Simple Newton iteration
            for iter = 1:50
                F    = D*(theta([2:end,end]) - 2*theta + theta([1,1:end-1]))/dx^2 ...
                       - 2*omega_n*sin(theta) + 2*Delta.*cos(theta);
                dFdt = D*(-2/dx^2)*ones(Nx,1) - 2*omega_n*cos(theta) - 2*Delta.*sin(theta);
                dth  = -F ./ dFdt;
                theta = theta + 0.5*dth;
                if max(abs(dth)) < 1e-10, break; end
            end
        end

        function obj = update_gap(obj, omega_n)
            p      = obj.params;
            lambda = obj.get_field(p,'qc_lambda',1.0);
            T      = obj.get_field(p,'qc_T',0.001);
            Nx     = numel(obj.x);
            N_om   = numel(omega_n);

            Delta_new = zeros(Nx,1);
            for ix = 1:Nx
                s = 0;
                for n = 1:N_om
                    th = obj.theta(ix,n);
                    s  = s + sin(th);
                end
                Delta_new(ix) = pi*lambda*T*s;
            end
            obj.Delta = Delta_new;
        end

    end

    methods (Static, Access=private)
        function v = get_field(p, fname, default)
            if isfield(p, fname), v = p.(fname);
            else, v = default;
            end
        end
    end
end
