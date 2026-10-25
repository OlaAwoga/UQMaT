classdef Spectral
    % UQMaT.Observables.Conventional.Spectral
    %   All conventional spectral observables in one class.
    %   Works for any system: normal, spinful, or superconducting (sc=true).
    %
    %   The Hamiltonian is diagonalised ONCE in the constructor.
    %   All subsequent method calls (dos, ldos, spin_texture, …) read
    %   from the cached eigenvalues E and eigenvectors V — no repeated eig().
    %
    %   obs = UQMaT.Observables.Conventional.Spectral(slv)
    %
    %   ---- Spectrum -----------------------------------------------
    %   [E, V]          = obs.spectrum()
    %
    %   ---- DOS ----------------------------------------------------
    %   [w, D]          = obs.dos()
    %   [w, D_up, D_dn] = obs.dos_spin()    % requires spin=true
    %   [w, D_e, D_h]   = obs.dos_ph()      % requires sc=true
    %
    %   ---- LDOS ---------------------------------------------------
    %   [w, L]          = obs.ldos()
    %   [w, L]          = obs.ldos(sites)
    %   [w, L_up, L_dn] = obs.ldos_spin(sites)
    %   [w, L_e, L_h]   = obs.ldos_ph(sites)
    %
    %   ---- Probability density ------------------------------------
    %   rho             = obs.prob_density()
    %   rho             = obs.prob_density(n_idx)
    %
    %   ---- Spin texture -------------------------------------------
    %   [w, Sx, Sy, Sz] = obs.spin_texture(sites)
    %
    %   ---- Pair amplitude -----------------------------------------
    %   F               = obs.pair_amplitude()  % sc=true only
    %
    %   ---- Energy grid --------------------------------------------
    %   All DOS/LDOS methods accept the energy grid in multiple forms:
    %
    %     obs.dos(-2:0.005:2)                    % positional vector
    %     obs.dos('EnergyGrid', -2:0.005:2)      % name-value
    %     obs.dos('range', [-2, 2])              % range → 500 pts auto
    %     obs.dos('range', [-2, 2], 'Npts', 300) % range + point count
    %     obs.dos('eta', 0.02)                   % broadening override
    %     obs.dos()                               % auto from spectrum
    %
    %   Priority: per-call argument > p.EnergyGrid > auto.

    %  Citation
    %  --------
    %  DOS via Lorentzian broadening: Haydock R. et al., J. Phys. C 5, 2845 (1972).
    %  Pair amplitude: de Gennes P.G., "Superconductivity of Metals and Alloys",
    %    Addison-Wesley (1966).
    %  Spin texture: Rashba E.I., Sov. Phys. Solid State 2, 1109 (1960).
    properties (SetAccess = private)
        solver              % UQMaT.Solver.Base
        N      (1,1) double % Ndim: physical sites
        p                   % UQMaT.Params
        E      (:,1) double % cached eigenvalues  (sorted)
        V      (:,:) double % cached eigenvectors (columns match E)
    end

    methods
        % ============================================================
        % Constructor — diagonalise once, cache E and V
        % ============================================================
        function obj = Spectral(solver)
            if ~isa(solver, 'UQMaT.Solver.Base')
                error('UQMaT:Conv:Spectral:badInput', ...
                    'Argument must be a UQMaT.Solver.Base subclass.');
            end
            obj.solver = solver;
            obj.N      = solver.sys.Ndim;
            obj.p      = solver.sys.params;

            % Single diagonalisation — all methods below read from cache
            [obj.E, obj.V] = solver.eigs();
        end

        % ============================================================
        % Spectrum
        % ============================================================
        function [E, V] = spectrum(obj)
            E = obj.E;
            V = obj.V;
        end

        % ============================================================
        % DOS — total density of states per physical site
        % ============================================================
        function [omega, D] = dos(obj, varargin)
            [omega, eta] = obj.grid(varargin{:});
            D = obj.broadened_sum(omega, eta, ones(size(obj.E))) / obj.N;
        end

        % ============================================================
        % DOS — spin-resolved  (requires spin=true)
        % ============================================================
        function [omega, D_up, D_dn] = dos_spin(obj, varargin)
            obj.needSpin('dos_spin');
            [omega, eta] = obj.grid(varargin{:});
            N   = obj.N;
            % Weight per state: total spin-up and spin-down electron weight
            w_up = sum(abs(obj.V(1:N,   :)).^2, 1)';   % Nstates×1
            w_dn = sum(abs(obj.V(N+1:2*N,:)).^2, 1)';
            D_up = obj.broadened_sum(omega, eta, w_up) / N;
            D_dn = obj.broadened_sum(omega, eta, w_dn) / N;
        end

        % ============================================================
        % DOS — particle/hole resolved  (requires sc=true)
        % ============================================================
        function [omega, D_e, D_h] = dos_ph(obj, varargin)
            obj.needSC('dos_ph');
            [omega, eta] = obj.grid(varargin{:});
            N   = obj.N;
            w_e = sum(abs(obj.V(1:2*N,     :)).^2, 1)';
            w_h = sum(abs(obj.V(2*N+1:4*N, :)).^2, 1)';
            D_e = obj.broadened_sum(omega, eta, w_e) / N;
            D_h = obj.broadened_sum(omega, eta, w_h) / N;
        end

        % ============================================================
        % LDOS — total local DOS  (electron block, summed over spin)
        % ============================================================
        function [omega, L] = ldos(obj, sites, varargin)
            if nargin < 2 || isempty(sites), sites = 1:obj.N; end
            [omega, eta] = obj.grid(varargin{:});
            N    = obj.N;
            Hdim = size(obj.V, 1);
            Ns   = numel(sites);
            Nw   = numel(omega);
            L    = zeros(Nw, Ns);
            % Block stride: how many rows separate spin-up from spin-down block.
            % For uniform systems: stride = N (normal/spinful/BdG all divide evenly).
            % For mixed heterostructures (e.g. SC substrate + normal chain),
            % Hdim is not a clean multiple of N, so we fall back to raw row weight.
            stride = round(Hdim / round(Hdim / N));  % nearest clean divisor
            use_spinsum = (Hdim >= 2*N) && (stride + max(sites) - 1 <= Hdim);
            for si = 1:Ns
                i = sites(si);
                if use_spinsum
                    w = abs(obj.V(i, :)).^2 + abs(obj.V(i+stride, :)).^2;
                else
                    w = abs(obj.V(i, :)).^2;
                end
                L(:,si) = obj.broadened_sum(omega, eta, w(:));
            end
        end

        % ============================================================
        % LDOS — spin-resolved  (requires spin=true)
        % ============================================================
        function [omega, L_up, L_dn] = ldos_spin(obj, sites, varargin)
            obj.needSpin('ldos_spin');
            if nargin < 2 || isempty(sites), sites = 1:obj.N; end
            [omega, eta] = obj.grid(varargin{:});
            N  = obj.N;
            Ns = numel(sites);
            Nw = numel(omega);
            L_up = zeros(Nw, Ns);
            L_dn = zeros(Nw, Ns);
            for si = 1:Ns
                i = sites(si);
                L_up(:,si) = obj.broadened_sum(omega, eta, abs(obj.V(i,   :)).^2');
                L_dn(:,si) = obj.broadened_sum(omega, eta, abs(obj.V(i+N, :)).^2');
            end
        end

        % ============================================================
        % LDOS — particle/hole resolved  (requires sc=true)
        % ============================================================
        function [omega, L_e, L_h] = ldos_ph(obj, sites, varargin)
            obj.needSC('ldos_ph');
            if nargin < 2 || isempty(sites), sites = 1:obj.N; end
            [omega, eta] = obj.grid(varargin{:});
            N    = obj.N;
            Hdim = size(obj.V, 1);
            Ns   = numel(sites);
            Nw   = numel(omega);
            L_e  = zeros(Nw, Ns);
            L_h  = zeros(Nw, Ns);

            % Determine block layout:
            % Uniform BdG (Hdim = 4N exactly): stride = N
            % Mixed heterostructure (Hdim ≠ 4N): stride = Hdim/4 rounded,
            %   but hole block may not exist for all sites — guard every access.
            if Hdim == 4*N
                stride = N;
            else
                stride = round(Hdim / 4);
            end

            for si = 1:Ns
                i  = sites(si);
                i2 = i + stride;      % spin-down electron
                i3 = i + 2*stride;    % spin-up hole
                i4 = i + 3*stride;    % spin-down hole

                % Electron weight (spin-up always present; spin-down if in range)
                if i2 <= Hdim
                    we = abs(obj.V(i,:)).^2 + abs(obj.V(i2,:)).^2;
                else
                    we = abs(obj.V(i,:)).^2;
                end

                % Hole weight (may not exist for mixed heterostructure sites)
                if i3 <= Hdim && i4 <= Hdim
                    wh = abs(obj.V(i3,:)).^2 + abs(obj.V(i4,:)).^2;
                elseif i3 <= Hdim
                    wh = abs(obj.V(i3,:)).^2;
                else
                    wh = zeros(1, size(obj.V,2));  % no hole block at this site
                end

                L_e(:,si) = obj.broadened_sum(omega, eta, we(:));
                L_h(:,si) = obj.broadened_sum(omega, eta, wh(:));
            end
        end

        % ============================================================
        % Probability density  |psi_n(i)|^2
        % ============================================================
        function rho = prob_density(obj, n_idx)
            N    = obj.N;
            Hdim = size(obj.V, 1);
            if nargin < 2 || isempty(n_idx)
                n_idx = 1:size(obj.V, 2);
            end
            Nst = numel(n_idx);
            rho = zeros(Nst, N);
            for k = 1:Nst
                n = n_idx(k);
                if Hdim == N
                    rho(k,:) = abs(obj.V(1:N, n)).^2;
                elseif Hdim == 2*N
                    rho(k,:) = abs(obj.V(1:N, n)).^2 + abs(obj.V(N+1:2*N, n)).^2;
                else  % 4N
                    rho(k,:) = abs(obj.V(1:N, n)).^2 + abs(obj.V(N+1:2*N, n)).^2;
                end
            end
        end

        % ============================================================
        % Spin texture  <S_alpha>(omega, site)  (requires spin=true)
        % ============================================================
        function [omega, Sx, Sy, Sz] = spin_texture(obj, sites, varargin)
            obj.needSpin('spin_texture');
            if nargin < 2 || isempty(sites), sites = 1:obj.N; end
            [omega, eta] = obj.grid(varargin{:});
            N   = obj.N;
            s   = UQMaT.SC.BdG.sigma();
            Nw  = numel(omega);
            Ns  = numel(sites);
            Sx  = zeros(Nw, Ns);
            Sy  = zeros(Nw, Ns);
            Sz  = zeros(Nw, Ns);

            % Pre-compute spinor at each site for all eigenstates (N_states)
            Nst = numel(obj.E);
            for si = 1:Ns
                i      = sites(si);
                % spinor: 2 × Nstates
                psi    = [obj.V(i,   :);
                          obj.V(i+N, :)];
                sx_n   = real(dot(psi, s.sx * psi, 1));   % 1 × Nstates
                sy_n   = real(dot(psi, s.sy * psi, 1));
                sz_n   = real(dot(psi, s.sz * psi, 1));
                Sx(:,si) = obj.broadened_sum(omega, eta, sx_n(:));
                Sy(:,si) = obj.broadened_sum(omega, eta, sy_n(:));
                Sz(:,si) = obj.broadened_sum(omega, eta, sz_n(:));
            end
        end

        % ============================================================
        % Pair amplitude  F(i) = sum_{E_n<0} u_n_up(i)*conj(v_n_dn(i))
        % (requires sc=true)
        % ============================================================
        function F = pair_amplitude(obj)
            obj.needSC('pair_amplitude');
            N    = obj.N;
            Hdim = size(obj.V, 1);
            F    = zeros(N, 1);
            occ  = find(obj.E < 0);
            % Block stride: standard BdG has Hdim=4N → stride=N.
            % Mixed heterostructure (e.g. SC substrate + normal chain):
            % Hdim ≠ 4N — use Hdim/4 as stride, clamped so indices stay valid.
            if Hdim == 4*N
                stride = N;
            else
                stride = floor(Hdim / 4);
            end
            % u_up: rows 1:N (electron spin-up)
            % v_dn: rows 3*stride+1 : 3*stride+N (hole spin-down)
            v_start = 3*stride + 1;
            v_end   = 3*stride + N;
            if v_end > Hdim
                % Hole block not fully available — use what exists
                v_end   = Hdim;
                Nv      = v_end - v_start + 1;
                for n = occ'
                    u_up = obj.V(1:Nv, n);
                    v_dn = obj.V(v_start:v_end, n);
                    F(1:Nv) = F(1:Nv) + u_up .* conj(v_dn);
                end
            else
                for n = occ'
                    u_up = obj.V(1:N,             n);
                    v_dn = obj.V(v_start:v_end,   n);
                    F    = F + u_up .* conj(v_dn);
                end
            end
        end

        % ============================================================
        function disp(obj)
            N = obj.N;
            p = obj.p;
            if p.sc
                mode = sprintf('SC BdG  (4×%d)', N);
            elseif p.spin
                mode = sprintf('spinful (2×%d)', N);
            else
                mode = sprintf('normal  (%d)',   N);
            end
            fprintf('  UQMaT.Observables.Conventional.Spectral\n');
            fprintf('  system=%s  %s\n', p.system, mode);
            fprintf('  E cached: %d states  [%.4g, %.4g]\n', ...
                numel(obj.E), obj.E(1), obj.E(end));
            fprintf('  Methods: spectrum, dos, dos_spin*, dos_ph**,\n');
            fprintf('           ldos, ldos_spin*, ldos_ph**, prob_density,\n');
            fprintf('           spin_texture*, pair_amplitude**\n');
            fprintf('  * requires spin=true  ** requires sc=true\n');
        end
    end

    % ---- Private helpers -------------------------------------------
    methods (Access = private)

        % Core broadening: D(omega) = sum_n w_n * Lorentzian(omega - E_n)
        % w  : Nstates×1  weight per eigenstate
        % returns Nomega×1 vector
        function D = broadened_sum(obj, omega, eta, w)
            % Vectorised: outer-product (Nomega × Nstates), then weighted sum
            % dE(i,n) = omega(i) - E(n)
            dE = omega(:) - obj.E(:)';          % Nomega × Nstates
            D  = (eta/pi) * (dE.^2 + eta^2).^(-1) * w(:);
        end

        function [omega, eta] = grid(obj, varargin)
            % Resolve the energy grid from multiple possible sources.
            %
            % Priority (first match wins):
            %   1. Positional vector:  obs.dos(-2:0.005:2)
            %   2. Name-value pair:    obs.dos('EnergyGrid', -2:0.005:2)
            %   3. Range shorthand:    obs.dos('range', [-2, 2])
            %   4. Range shorthand:    obs.dos('range', [-2, 2], 'Npts', 300)
            %   5. 'eta' override:     obs.dos('eta', 0.02)
            %   6. p.EnergyGrid        set in Params
            %   7. auto                from spectrum extent ± padding (500 pts)
            p   = obj.p;
            eta = p.eta;
            omega = [];

            % Parse all name-value and positional options
            k = 1;
            while k <= numel(varargin)
                arg = varargin{k};
                if isnumeric(arg) && isvector(arg) && numel(arg) > 2
                    % Positional energy vector: obs.dos(-2:0.005:2)
                    omega = arg(:);
                    k = k + 1;
                elseif ischar(arg) || isstring(arg)
                    key = lower(strtrim(char(arg)));
                    val = varargin{k+1};
                    switch key
                        case 'energygrid'
                            omega = val(:);
                        case 'range'
                            % 'range', [Emin, Emax]  — Npts from next arg or 500
                            Emin = val(1);  Emax = val(2);
                            Npts = 500;
                            if k+2 <= numel(varargin) && ischar(varargin{k+2}) && ...
                                    strcmpi(varargin{k+2},'Npts')
                                Npts = varargin{k+3};
                                k    = k + 2;
                            end
                            omega = linspace(Emin, Emax, Npts)';
                        case 'npts'
                            % 'Npts', 300  — set number of auto points
                            E   = obj.E;
                            pad = max(0.2*(E(end)-E(1)), 4*eta);
                            omega = linspace(E(1)-pad, E(end)+pad, val)';
                        case 'eta'
                            eta = val;
                        case 'emin'
                            % 'emin', -2, 'emax', 2
                            Emin = val;
                            emax_idx = find(cellfun(@(x) ischar(x) && ...
                                strcmpi(x,'emax'), varargin(k+2:end)), 1);
                            if ~isempty(emax_idx)
                                Emax = varargin{k+2+emax_idx};
                            else
                                Emax = obj.E(end) + 0.1;
                            end
                            omega = linspace(Emin, Emax, 500)';
                    end
                    k = k + 2;
                else
                    k = k + 1;
                end
            end

            % Fallback: Params or auto
            if isempty(omega)
                if ~isempty(p.EnergyGrid)
                    omega = p.EnergyGrid(:);
                else
                    E   = obj.E;
                    pad = max(0.2*(E(end)-E(1)), 4*eta);
                    omega = linspace(E(1)-pad, E(end)+pad, 500)';
                end
            end
        end

        function needSpin(obj, mname)
            if isempty(obj.p), return; end
            if obj.p.spin, return; end
            % Heterostructure: check if any side has spin-active fields
            p = obj.p;
            if ~isempty(p.hetero) && isstruct(p.hetero)
                h = p.hetero;
                sides = {};
                if isfield(h,'substrate'), sides{end+1} = h.substrate; end
                if isfield(h,'L'),         sides{end+1} = h.L; end
                if isfield(h,'R'),         sides{end+1} = h.R; end
                if isfield(h,'chain'),     sides{end+1} = h.chain; end
                for k = 1:numel(sides)
                    s = sides{k};
                    if isfield(s,'B') && any(s.B ~= 0), return; end
                    if isfield(s,'lambda_R') && s.lambda_R ~= 0, return; end
                    if isfield(s,'SC') && ~isempty(fieldnames(s.SC)), return; end
                end
            end
            error('UQMaT:Conv:Spectral:noSpin', ...
                '%s requires spin=true (set B, lambda_R, or sc=true).', mname);
        end

        function needSC(obj, mname)
            if isempty(obj.p), return; end  % no Params — assume SC (fromEV path)
            if obj.p.sc, return; end        % explicit sc=true — pass
            % Heterostructure: check if any side has SC pairing
            p = obj.p;
            if ~isempty(p.hetero) && isstruct(p.hetero)
                h = p.hetero;
                sides = {};
                if isfield(h,'substrate'), sides{end+1} = h.substrate; end
                if isfield(h,'L'),         sides{end+1} = h.L; end
                if isfield(h,'R'),         sides{end+1} = h.R; end
                if isfield(h,'chain'),     sides{end+1} = h.chain; end
                for k = 1:numel(sides)
                    s = sides{k};
                    if isfield(s,'SC') && ~isempty(fieldnames(s.SC))
                        return   % heterostructure has SC component — pass
                    end
                    if isfield(s,'Delta_s') && any(s.Delta_s ~= 0)
                        return
                    end
                    if isfield(s,'Delta_p') && any(s.Delta_p ~= 0)
                        return
                    end
                end
            end
            error('UQMaT:Conv:Spectral:noSC', ...
                '%s requires sc=true (set Delta_s, dvec, or Delta_p nonzero).', mname);
        end
    end

    methods (Static)
        function obj = fromEV(ev)
            % Construct Spectral directly from a UQMaT.EV cache.
            % No solver required — reads E and V from the cache.
            %   obs = UQMaT.Observables.Conventional.Spectral.fromEV(ev)
            %
            % This is called automatically by ev.spectral().
            % The solver field is left empty; all methods read from E,V.
            obj        = UQMaT.Observables.Conventional.Spectral.make_empty_();
            obj.N      = ev.N;
            obj.p      = ev.p;
            obj.E      = ev.E;
            obj.V      = ev.V;
        end
    end

    methods (Static, Access = private)
        function obj = make_empty_()
            % Internal: create an uninitialised Spectral object.
            % Bypasses the solver check by going through a minimal
            % constructor path — solver is set to [] (never used after fromEV).
            p   = UQMaT.Params('system','chain','N',2,'t',1,'mu',0);
            sys = UQMaT.System(p);
            slv = UQMaT.Solver.Diagonalization.Direct(sys);
            obj = UQMaT.Observables.Conventional.Spectral(slv);
            % Fields will be overwritten by fromEV — the tiny 2-site
            % system is diagonalised in microseconds and immediately
            % discarded. The constructor always runs once here.
        end
    end
end

        % ============================================================
        % Self-consistent BdG gap equation
        % ============================================================
        function [Delta_conv, F_conv, history] = selfconsistent(obj, lambda, varargin)
            % Iterate the gap equation  Δ(i) = λ · F(i)  to self-consistency.
            %
            %   F(i) = Σ_{E_n < 0}  u_{n↑}(i) · conj(v_{n↓}(i))
            %   Δ(i) = λ · F(i)
            %
            %   The Hamiltonian is rebuilt each iteration with the updated
            %   Delta_s_vec. Convergence is declared when
            %   max|Δ_new - Δ_old| < tol.
            %
            %   [Delta, F, history] = obs.selfconsistent(lambda)
            %   [Delta, F, history] = obs.selfconsistent(lambda, 'tol', 1e-6)
            %   [Delta, F, history] = obs.selfconsistent(lambda, 'maxiter', 200)
            %   [Delta, F, history] = obs.selfconsistent(lambda, 'mix', 0.5)
            %     mix : linear mixing parameter  Δ_new = mix*Δ_iter + (1-mix)*Δ_old
            %           (0 < mix ≤ 1, default 0.5; smaller = more stable but slower)
            %   [Delta, F, history] = obs.selfconsistent(lambda, 'verbose', true)
            %
            %   Returns:
            %     Delta_conv : N×1  converged site-resolved gap Δ(i)
            %     F_conv     : N×1  converged pair amplitude F(i)
            %     history    : struct with fields:
            %                    .Delta  : Niter×N  gap at each iteration
            %                    .error  : Niter×1  max|ΔΔ| at each step
            %                    .converged : logical

            obj.needSC('selfconsistent');

            % Options
            opt = struct('tol',1e-6,'maxiter',200,'mix',0.5,'verbose',false);
            for k = 1:2:numel(varargin)
                opt.(varargin{k}) = varargin{k+1};
            end

            p   = obj.p;
            sys = obj.solver.sys;
            N   = obj.N;

            % Initial Delta from current Params
            Delta = p.Delta_s_vec(1:N);

            history.Delta     = zeros(opt.maxiter, N);
            history.error     = zeros(opt.maxiter, 1);
            history.converged = false;

            if opt.verbose
                fprintf('Self-consistent BdG: λ=%.4g  tol=%.2e  maxiter=%d\n', ...
                    lambda, opt.tol, opt.maxiter);
            end

            for iter = 1:opt.maxiter
                % ---- Rebuild Params with current Delta ----
                p_new = UQMaT.Params('system', p.system, ...
                    'N',         p.N,          ...
                    'bc',        p.bc,         ...
                    't',         p.t,          ...
                    'mu',        p.mu,         ...
                    't2',        p.t2,         ...
                    'onsite',    p.onsite,     ...
                    'B',         p.B,          ...
                    'lambda_R',  p.lambda_R,   ...
                    'sc',        true,         ...
                    'Delta_s',   Delta,        ...
                    'W',         p.W,          ...
                    'seed',      p.seed,       ...
                    'eta',       p.eta,        ...
                    'EnergyGrid', p.EnergyGrid);

                sys_new = UQMaT.System(p_new);
                slv_new = UQMaT.Solver.Diagonalization.Direct(sys_new);
                obs_new = UQMaT.Observables.Conventional.Spectral(slv_new);

                % ---- Compute pair amplitude F(i) ----
                F_new = obs_new.pair_amplitude();

                % ---- Update gap ----
                Delta_new = lambda * F_new;

                % ---- Linear mixing ----
                Delta_mix = opt.mix * Delta_new + (1 - opt.mix) * Delta;

                % ---- Convergence check ----
                err = max(abs(Delta_mix - Delta));
                history.Delta(iter,:) = Delta_mix';
                history.error(iter)   = err;

                if opt.verbose
                    fprintf('  iter %3d: max|ΔΔ| = %.3e  mean|Δ| = %.4f\n', ...
                        iter, err, mean(abs(Delta_mix)));
                end

                Delta = Delta_mix;

                if err < opt.tol
                    history.converged = true;
                    history.Delta     = history.Delta(1:iter,:);
                    history.error     = history.error(1:iter);
                    if opt.verbose
                        fprintf('  Converged at iteration %d.\n', iter);
                    end
                    break
                end
            end

            if ~history.converged
                warning('UQMaT:Spectral:noConvergence', ...
                    'Self-consistent gap did not converge after %d iterations (tol=%.2e, final err=%.2e).', ...
                    opt.maxiter, opt.tol, err);
                history.Delta = history.Delta(1:opt.maxiter,:);
                history.error = history.error(1:opt.maxiter);
            end

            Delta_conv = Delta;
            F_conv     = F_new;

            % ---- Update internal cache with converged solution ----
            % Rebuild one final time so obs can be used for further analysis
            p_fin = UQMaT.Params('system', p.system, ...
                'N',        p.N,          'bc',      p.bc,    ...
                't',        p.t,          'mu',      p.mu,    ...
                't2',       p.t2,         'onsite',  p.onsite, ...
                'B',        p.B,          'lambda_R',p.lambda_R, ...
                'sc',       true,         'Delta_s', Delta_conv, ...
                'W',        p.W, ...
                'seed',     p.seed,       'eta',     p.eta, ...
                'EnergyGrid', p.EnergyGrid);
            sys_fin     = UQMaT.System(p_fin);
            slv_fin     = UQMaT.Solver.Diagonalization.Direct(sys_fin);
            [E_fin, V_fin] = slv_fin.eigs();
            obj.E = E_fin;
            obj.V = V_fin;
            obj.p = p_fin;
        end
