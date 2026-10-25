classdef EV < handle
    % UQMaT.EV  Single diagonalisation result shared across all observables.
    %
    %   Solve once, pass to any observable — no repeated solve calls.
    %
    %   ================================================================
    %   Construction — choose your solver
    %   ================================================================
    %
    %   Default (dense diagonalisation, all states):
    %     ev = UQMaT.EV(sys)
    %
    %   Explicit solver name (string):
    %     ev = UQMaT.EV(sys, 'direct')              % dense eig(), all states
    %     ev = UQMaT.EV(sys, 'lanczos')             % sparse eigs(), k=20 near lm
    %     ev = UQMaT.EV(sys, 'lanczos', 'k',10, 'sigma',0)  % 10 states near E=0
    %     ev = UQMaT.EV(sys, 'resolvent')           % G(ω) — DOS via inversion
    %     ev = UQMaT.EV(sys, 'resolvent', 'eta',0.02)
    %     ev = UQMaT.EV(sys, 'chebyshev')           % KPM — large N, no eig
    %     ev = UQMaT.EV(sys, 'chebyshev', 'M',512)
    %
    %   Pass a pre-built solver object:
    %     slv = UQMaT.Solver.Diagonalization.Lanczos(sys, 'k',10, 'sigma',0);
    %     ev  = UQMaT.EV(slv)
    %
    %     slv = UQMaT.Solver.Greens.Resolvent(sys, 'eta',0.01);
    %     ev  = UQMaT.EV(slv)
    %
    %     slv = UQMaT.Solver.Spectral.Chebyshev(sys, 'M',512);
    %     ev  = UQMaT.EV(slv)
    %
    %   Supply your own V and E directly:
    %     ev = UQMaT.EV(V, E)           % Hdim×Nstates, Nstates×1
    %     ev = UQMaT.EV(V, E, N)        % with explicit site count
    %     ev = UQMaT.EV(V, E, N, p)     % with Params
    %
    %   ================================================================
    %   Solver reference
    %   ================================================================
    %
    %   Solver              Class                             Best for
    %   ─────────────────────────────────────────────────────────────────
    %   'direct'     Solver.Diagonalization.Direct    N ≲ 2000, all states
    %   'lanczos'    Solver.Diagonalization.Lanczos   Large N, few states
    %   'resolvent'  Solver.Greens.Resolvent           G(ω) DOS via inversion
    %   'chebyshev'  Solver.Spectral.Chebyshev         Very large N, KPM DOS
    %
    %   Solver options (all via name-value after solver name):
    %
    %   Lanczos:
    %     'k'      number of eigenvalues          (default 20)
    %     'sigma'  shift — eigenvalues near sigma  (default 'lm' = largest mag.)
    %              pass 0 for states near E=0
    %     'tol'    convergence tolerance           (default 0)
    %     'maxit'  max iterations                  (default 300)
    %
    %   Resolvent:
    %     'eta'    Lorentzian broadening           (default from p.eta)
    %     'EnergyGrid' or positional vector
    %
    %   Chebyshev:
    %     'M'      number of Chebyshev moments     (default 256)
    %     'Nrand'  random vectors for stochastic trace (default 10)
    %     'kernel' 'jackson'|'lorentz'             (default 'jackson')
    %
    %   ================================================================
    %   Properties
    %   ================================================================
    %     ev.E      : Nstates×1  sorted eigenvalues
    %     ev.V      : Hdim×Nstates  eigenvectors (columns sorted by E)
    %     ev.N      : number of physical sites (Ndim)
    %     ev.p      : UQMaT.Params
    %     ev.solver : the solver object used
    %
    %   ================================================================
    %   Observables — all read from cached ev.E and ev.V
    %   ================================================================
    %     obs = ev.spectral()              % Conventional.Spectral
    %     inv = ev.realspace_invariants()  % Topological.RealSpace.Invariants
    %     mp  = ev.mp()                    % RealSpace.MajoranaPolarization
    %
    %   ================================================================
    %   Convenience pass-throughs (energy range in all forms)
    %   ================================================================
    %     [E, V]          = ev.spectrum()
    %     idx             = ev.edge_idx()
    %
    %     [w, D]          = ev.dos()
    %     [w, D]          = ev.dos(-2:0.005:2)
    %     [w, D]          = ev.dos('range',[-2,2])
    %     [w, D]          = ev.dos('range',[-2,2],'Npts',300)
    %     [w, D]          = ev.dos('eta',0.02)
    %     [w, D_up, D_dn] = ev.dos_spin()
    %     [w, D_e, D_h]   = ev.dos_ph()
    %     [w, L]          = ev.ldos(sites)
    %     [w, L_up, L_dn] = ev.ldos_spin(sites)
    %     [w, L_e, L_h]   = ev.ldos_ph(sites)
    %     [w, Sx, Sy, Sz] = ev.spin_texture(sites)
    %     F               = ev.pair_amplitude()
    %     rho             = ev.prob_density(n_idx)
    %
    %  Citation
    %  --------
    %  LAPACK eig: Anderson E. et al., LAPACK Users' Guide, SIAM (1999).
    %  Lanczos/ARPACK: Lehoucq R.B. et al., ARPACK Users' Guide, SIAM (1998).
    %  Green's function resolvent: Economou E.N., Green's Functions in
    %    Quantum Physics, Springer (2006).
    %  KPM: Weisse A. et al., Rev. Mod. Phys. 78, 275 (2006).

    properties (SetAccess = private)
        E      (:,1) double     % sorted eigenvalues
        V      (:,:) double     % eigenvectors, columns sorted by E
        N      (1,1) double     % number of physical sites (Ndim)
        Hdim   (1,1) double     % Hamiltonian dimension
        p                       % UQMaT.Params (may be empty)
        solver                  % the solver object used
        solver_name (1,:) char  % 'direct'|'lanczos'|'resolvent'|'chebyshev'|'external'
        % Lazily constructed observable handles (handle class — mutations persist)
        obs_spectral_ = []
        obs_mp_       = []
    end

    methods

        function obj = EV(varargin)
            if nargin == 0
                error('UQMaT:EV:noArgs', 'EV requires at least one argument.');
            end

            first = varargin{1};

            % ── From raw eigenvectors ──────────────────────────────────────
            % ── From raw eigenvectors:  EV(V, E)  or  EV(V, E, N)  or  EV(V, E, N, p) ──
            if isnumeric(first)
                if nargin < 2 || ~isnumeric(varargin{2})
                    error('UQMaT:EV:badArgs', ...
                        'EV(V, E): second argument must be a numeric eigenvalue vector.');
                end
                V_in = first;
                E_in = varargin{2};
                [obj.E, idx] = sort(real(E_in(:)));
                obj.V    = V_in(:, idx);
                obj.Hdim = size(obj.V, 1);
                if nargin >= 3 && ~isempty(varargin{3})
                    obj.N = varargin{3};
                else
                    obj.N = obj.Hdim;
                end
                if nargin >= 4 && ~isempty(varargin{4})
                    obj.p = varargin{4};
                else
                    obj.p = [];
                end
                obj.solver_name = 'external';
                return
            end

            % ── From a pre-built solver object:  EV(slv) ──────────────────
            if isa(first, 'UQMaT.Solver.Base')
                obj = obj.run_solver_(first, class(first));
                return
            end

            % ── From System:  EV(sys)  or  EV(sys, 'solver', opts...) ────
            if ~isa(first, 'UQMaT.System')
                error('UQMaT:EV:badArgs', ...
                    ['First argument must be UQMaT.System, UQMaT.Solver.Base, ' ...
                     'or a numeric matrix V.\n' ...
                     'Usage:\n' ...
                     '  EV(sys)                         — dense (default)\n' ...
                     '  EV(sys, ''direct'')               — dense eig()\n' ...
                     '  EV(sys, ''lanczos'', ''k'',10, ''sigma'',0)\n' ...
                     '  EV(sys, ''resolvent'', ''eta'',0.02)\n' ...
                     '  EV(sys, ''chebyshev'', ''M'',512)\n' ...
                     '  EV(slv)                         — pre-built solver\n' ...
                     '  EV(V, E)                        — raw eigenvectors\n' ...
                     '  EV(V, E, N, p)                  — with site count and Params']);
            end

            sys = first;

            % No second argument → default dense solver
            if nargin < 2 || isempty(varargin{2})
                slv = UQMaT.Solver.Diagonalization.Direct(sys);
                obj = obj.run_solver_(slv, 'direct');
                return
            end

            % Second argument is a solver name string
            if ischar(varargin{2}) || isstring(varargin{2})
                sname = lower(strtrim(char(varargin{2})));
                opts  = varargin(3:end);
                slv   = obj.build_solver_(sys, sname, opts);
                obj   = obj.run_solver_(slv, sname);
                return
            end

            % Second argument is unrecognised
            error('UQMaT:EV:badSolverArg', ...
                ['Second argument must be a solver name string:\n' ...
                 '  ''direct'' | ''lanczos'' | ''resolvent'' | ''chebyshev''\n' ...
                 'Got: %s'], class(varargin{2}));

            error('UQMaT:EV:badSolverArg', ...
                ['Second argument must be a solver name string:\n' ...
                 '  ''direct'' | ''lanczos'' | ''resolvent'' | ''chebyshev''\n' ...
                 'or omit it for the default dense solver.']);
        end

        % ============================================================
        % In-gap state indices
        % ============================================================
        function idx = edge_idx(obj, tol)
            % Indices of in-gap (near-zero energy) eigenstates.
            %   idx = ev.edge_idx()        % auto threshold
            %   idx = ev.edge_idx(0.05)    % explicit threshold
            if isempty(obj.E)
                idx = []; return
            end
            if nargin < 2 || isempty(tol)
                pos = obj.E(obj.E >  1e-10);
                neg = obj.E(obj.E < -1e-10);
                if ~isempty(pos) && ~isempty(neg)
                    tol = min(pos(1), abs(neg(end))) * 0.5;
                else
                    tol = 0.05 * (obj.E(end) - obj.E(1));
                end
            end
            idx = find(abs(obj.E) < tol);
        end

        % ============================================================
        % Observable handles (lazily built, cached by handle semantics)
        % ============================================================
        function obs = spectral(obj)
            if isempty(obj.obs_spectral_)
                obj.obs_spectral_ = ...
                    UQMaT.Observables.Conventional.Spectral.fromEV(obj);
            end
            obs = obj.obs_spectral_;
        end

        function mp = mp(obj)
            if isempty(obj.obs_mp_)
                obj.obs_mp_ = ...
                    UQMaT.Observables.Topological.RealSpace.MajoranaPolarization(...
                    obj.V, obj.N, obj.E);
            end
            mp = obj.obs_mp_;
        end

        function inv = realspace_invariants(obj)
            inv = UQMaT.Observables.Topological.RealSpace.Invariants.fromEV(obj);
        end

        % ============================================================
        % DOS / LDOS pass-throughs
        % All energy arguments forwarded to Spectral.grid():
        %   ev.dos(-2:0.005:2)
        %   ev.dos('range',[-2,2])
        %   ev.dos('range',[-2,2],'Npts',300)
        %   ev.dos('eta',0.02)
        %   ev.dos()
        % ============================================================
        function [w, D] = dos(obj, varargin)
            [w, D] = obj.spectral().dos(varargin{:});
        end

        function [w, D_up, D_dn] = dos_spin(obj, varargin)
            [w, D_up, D_dn] = obj.spectral().dos_spin(varargin{:});
        end

        function [w, D_e, D_h] = dos_ph(obj, varargin)
            [w, D_e, D_h] = obj.spectral().dos_ph(varargin{:});
        end

        function [w, L] = ldos(obj, sites, varargin)
            [w, L] = obj.spectral().ldos(sites, varargin{:});
        end

        function [w, L_up, L_dn] = ldos_spin(obj, sites, varargin)
            [w, L_up, L_dn] = obj.spectral().ldos_spin(sites, varargin{:});
        end

        function [w, L_e, L_h] = ldos_ph(obj, sites, varargin)
            [w, L_e, L_h] = obj.spectral().ldos_ph(sites, varargin{:});
        end

        function [w, Sx, Sy, Sz] = spin_texture(obj, sites, varargin)
            [w, Sx, Sy, Sz] = obj.spectral().spin_texture(sites, varargin{:});
        end

        function F = pair_amplitude(obj)
            F = obj.spectral().pair_amplitude();
        end

        function rho = prob_density(obj, n_idx)
            rho = obj.spectral().prob_density(n_idx);
        end

        function [E, V] = spectrum(obj)
            E = obj.E;  V = obj.V;
        end

        % ============================================================
        function disp(obj)
            fprintf('\n  UQMaT.EV\n');
            fprintf('  Solver  : %s\n', obj.solver_name);
            if ~isempty(obj.p)
                fprintf('  System  : %s   N=%d   Hdim=%d\n', ...
                    obj.p.system, obj.N, obj.Hdim);
            else
                fprintf('  N=%d   Hdim=%d\n', obj.N, obj.Hdim);
            end
            fprintf('  Nstates : %d\n', numel(obj.E));
            if ~isempty(obj.E)
                fprintf('  E range : [%.4g, %.4g]\n', obj.E(1), obj.E(end));
                idx = obj.edge_idx();
                if ~isempty(idx)
                    fprintf('  In-gap  : %d states  E = [', numel(idx));
                    fprintf(' %.3g', obj.E(idx)); fprintf(' ]\n');
                else
                    fprintf('  In-gap  : none\n');
                end
            else
                fprintf('  E : (not cached — use dos/ldos directly on solver)\n');
            end
            fprintf('\n');
        end

    end   % public methods

    % ================================================================
    % Private helpers
    % ================================================================
    methods (Access = private)

        function slv = build_solver_(~, sys, sname, opts)
            % Build a solver object from a name string and option pairs.
            switch sname
                case 'direct'
                    slv = UQMaT.Solver.Diagonalization.Direct(sys);

                case 'lanczos'
                    slv = UQMaT.Solver.Diagonalization.Lanczos(sys);
                    for i = 1:2:numel(opts)
                        slv.(opts{i}) = opts{i+1};
                    end

                case 'resolvent'
                    slv = UQMaT.Solver.Greens.Resolvent(sys);
                    for i = 1:2:numel(opts)
                        slv.(opts{i}) = opts{i+1};
                    end

                case 'chebyshev'
                    slv = UQMaT.Solver.Spectral.Chebyshev(sys);
                    for i = 1:2:numel(opts)
                        slv.(opts{i}) = opts{i+1};
                    end

                otherwise
                    error('UQMaT:EV:unknownSolver', ...
                        ['Unknown solver: "%s"\n' ...
                         'Valid names: ''direct'' | ''lanczos'' | ' ...
                         '''resolvent'' | ''chebyshev'''], sname);
            end
        end

        function obj = run_solver_(obj, slv, sname)
            % Run solver, cache E/V, set metadata.
            obj.solver      = slv;
            obj.solver_name = sname;
            if isprop(slv,'sys') || isfield(slv,'sys')
                sys      = slv.sys;
                obj.N    = sys.Ndim;
                obj.p    = sys.params;
                obj.Hdim = size(sys.H, 1);
            end
            % All solvers expose eigs() returning [E, V]
            [obj.E, obj.V] = slv.eigs();
        end

    end   % private methods

end
