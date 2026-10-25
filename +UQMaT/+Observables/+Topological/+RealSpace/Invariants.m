classdef Invariants
    % UQMaT.Observables.Topological.RealSpace.Invariants
    %   Real-space (OBC) signatures of topological phases.
    %
    %   Construction:
    %     inv = UQMaT.Observables.Topological.RealSpace.Invariants(sys)
    %     inv = ev.realspace_invariants()   % from EV — no re-diagonalisation
    %
    %   Methods:
    %     [E_edge, idx] = inv.edge_states(tol)
    %     [rho, x]      = inv.edge_density(tol)
    %     [ipr, E]      = inv.ipr()
    %     F             = inv.pair_amplitude()
    %     [F, x]        = inv.pair_amplitude_profile()
    %     mp            = inv.majorana_polarization()
    %     disp(inv)
    %
    %  Citation
    %  --------
    %  Edge states / bulk-edge correspondence:
    %    Hatsugai Y., PRL 71, 3697 (1993).
    %  IPR (inverse participation ratio):
    %    Bell R.J., Dean P., Disc. Faraday Soc. 50, 55 (1970).
    %  Pair amplitude at edges:
    %    Zhu J.-X., BdG Method, Springer (2016). [Chapter 4]

    properties (SetAccess = private)
        obs  UQMaT.Observables.Conventional.Spectral
        % Cached from EV or System — used when sys is unavailable
        N_   (1,1) double = 0
        p_                      % UQMaT.Params (may be empty)
        E_   (:,1) double = []  % eigenvalues cache
    end

    methods (Static)
        function obj = fromEV(ev)
            % Construct from a UQMaT.EV cache — no re-diagonalisation.
            %   inv = UQMaT.Observables.Topological.RealSpace.Invariants.fromEV(ev)
            obj      = UQMaT.Observables.Topological.RealSpace.Invariants.make_empty_();
            obj.obs  = UQMaT.Observables.Conventional.Spectral.fromEV(ev);
            obj.N_   = ev.N;
            obj.p_   = ev.p;
            obj.E_   = ev.E;
        end
    end

    methods (Static, Access = private)
        function obj = make_empty_()
            p   = UQMaT.Params('system','chain','N',2,'t',1,'mu',0);
            sys = UQMaT.System(p);
            obj = UQMaT.Observables.Topological.RealSpace.Invariants(sys);
        end
    end

    methods
        function obj = Invariants(sys)
            if ~isa(sys, 'UQMaT.System')
                error('UQMaT:Topo:RS:badInput','Argument must be a UQMaT.System.');
            end
            slv     = UQMaT.Solver.Diagonalization.Direct(sys);
            obj.obs = UQMaT.Observables.Conventional.Spectral(slv);
            obj.N_  = sys.Ndim;
            obj.p_  = sys.params;
            obj.E_  = sys.E;
        end

        % ============================================================
        % Edge / in-gap states
        % ============================================================
        function [E_edge, idx] = edge_states(obj, tol)
            E = obj.E_;
            if nargin < 2 || isempty(tol)
                pos = E(E >  1e-10);
                neg = E(E < -1e-10);
                if ~isempty(pos) && ~isempty(neg)
                    tol = min(pos(1), abs(neg(end))) * 0.5;
                else
                    tol = 0.05 * (E(end) - E(1));
                end
            end
            idx    = find(abs(E) < tol);
            E_edge = E(idx);
        end

        % ============================================================
        % Probability density of edge / in-gap states
        % ============================================================
        function [rho, x] = edge_density(obj, tol)
            if nargin < 2 || isempty(tol)
                [~, idx] = obj.edge_states();
            else
                [~, idx] = obj.edge_states(tol);
            end
            if isempty(idx)
                warning('UQMaT:Topo:RS:noEdge','No in-gap states found. Try increasing tol.');
                rho = []; x = []; return
            end
            rho = obj.obs.prob_density(idx);   % Nedge × N
            x   = (1:obj.N_)';
        end

        % ============================================================
        % Inverse Participation Ratio
        % ============================================================
        function [ipr, E] = ipr(obj)
            [E, V] = obj.obs.spectrum();
            N    = obj.N_;
            Hdim = size(V, 1);
            Nst  = size(V, 2);
            ipr  = zeros(Nst, 1);
            for n = 1:Nst
                if Hdim == N
                    psi = abs(V(1:N, n)).^2;
                elseif Hdim >= 2*N
                    psi = abs(V(1:N, n)).^2 + abs(V(N+1:2*N, n)).^2;
                end
                ipr(n) = sum(psi.^2);
            end
        end

        % ============================================================
        % Pair amplitude
        % ============================================================
        function F = pair_amplitude(obj)
            F = obj.obs.pair_amplitude();
        end

        function [F, x] = pair_amplitude_profile(obj)
            F = obj.obs.pair_amplitude();
            x = (1:obj.N_)';
        end

        % ============================================================
        % Majorana Polarization
        % ============================================================
        function mp = majorana_polarization(obj)
            % Return a MajoranaPolarization object.
            %   mp = inv.majorana_polarization()
            [E, V] = obj.obs.spectrum();
            mp = UQMaT.Observables.Topological.RealSpace.MajoranaPolarization(...
                V, obj.N_, E);
        end

        % ============================================================
        function disp(obj)
            p = obj.p_;
            E = obj.E_;
            fprintf('\n  UQMaT.Observables.Topological.RealSpace.Invariants\n');
            if ~isempty(p)
                fprintf('  system=%s  N=%d  bc=%s\n\n', p.system, obj.N_, upper(p.bc));
            else
                fprintf('  N=%d\n\n', obj.N_);
            end

            [E_edge, idx] = obj.edge_states();
            if ~isempty(E_edge)
                fprintf('  In-gap states: %d found\n', numel(E_edge));
                fprintf('  Energies: '); fprintf('%.4g  ', E_edge); fprintf('\n');
            else
                fprintf('  In-gap states: none (fully gapped or trivial)\n');
            end

            if ~isempty(idx)
                [ipr_all, ~] = obj.ipr();
                ipr_edge = ipr_all(idx);
                fprintf('  IPR of edge states: ');
                fprintf('%.4g  ', ipr_edge);
                fprintf('  (extended=%.4g, localised=1)\n', 1/obj.N_);
            end

            % Majorana Polarization
            is_bdg = ~isempty(p) && (p.sc || ismember(p.system,{'kitaev','hnkitaev'}));
            if is_bdg && ~isempty(idx)
                try
                    mp = obj.majorana_polarization();
                    [~, E_mp, idx_mp] = mp.local_edge();
                    if ~isempty(idx_mp)
                        [P, M]           = mp.total(idx_mp);
                        [MP_L,MP_R,MP_LR] = mp.left_right(idx_mp);
                        fprintf('  Majorana Polarization:\n');
                        for k = 1:numel(idx_mp)
                            fprintf('    state %d (E=%.4g): P=%.4f  M=%.4f  MP_L=%.4f  MP_R=%.4f  MP_LR=%.4f\n', ...
                                idx_mp(k), E_mp(k), P(k), M(k), MP_L(k), MP_R(k), MP_LR(k));
                        end
                    end
                catch
                end
            end

            % Pair amplitude (SC only, not Kitaev)
            if is_bdg && ~isempty(p) && ~strcmp(p.system,'kitaev')
                try
                    F = obj.pair_amplitude();
                    fprintf('  |F| max=%.4g  mean=%.4g\n', max(abs(F)), mean(abs(F)));
                catch
                end
            end

            pos = E(E >  1e-10);
            neg = E(E < -1e-10);
            if ~isempty(pos) && ~isempty(neg)
                fprintf('  Quasiparticle gap: %.6g\n', pos(1) - neg(end));
            else
                fprintf('  Quasiparticle gap: gapless\n');
            end
            fprintf('\n');
        end
    end
end
