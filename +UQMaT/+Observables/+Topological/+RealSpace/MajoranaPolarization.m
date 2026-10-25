classdef MajoranaPolarization
    % UQMaT.Observables.Topological.RealSpace.MajoranaPolarization
    %   Majorana Polarization and Charge of BdG eigenstates.
    %
    %   BdG eigenvector (spinful, 4N basis):
    %     ψ_n = ( u_{n↑}(i), u_{n↓}(i), v_{n↑}(i), v_{n↓}(i) )ᵀ
    %   BdG eigenvector (spinless, 2N basis):
    %     ψ_n = ( u_n(i), v_n(i) )ᵀ
    %
    %   ================================================================
    %   CHARGE  C_n(i)
    %   ================================================================
    %   Local charge (electron-hole imbalance) at site i:
    %
    %     C_n(i) = |u_n(i)|² - |v_n(i)|²   ∈ [-1, +1]
    %
    %   where  |u_n(i)|² = Σ_σ |u_{nσ}(i)|²  (electron weight, spinful)
    %          |v_n(i)|² = Σ_σ |v_{nσ}(i)|²  (hole weight, spinful)
    %
    %   C_n(i) = +1 : purely electronic state
    %   C_n(i) = -1 : purely hole state
    %   C_n(i) =  0 : equal electron/hole weight (Majorana condition)
    %
    %   ================================================================
    %   LOCAL MAJORANA POLARIZATION  MP_n(i)
    %   ================================================================
    %   Spinless BdG (2N basis):
    %     MP_n(i) = u_n(i)·v_n*(i) + u_n*(i)·v_n(i)
    %             = 2 Re[ u_n(i)·v_n*(i) ]
    %
    %   Spinful BdG (4N basis):
    %     MP_n(i) = 2 Σ_σ Re[ u_{nσ}(i)·v_{nσ}*(i) ]
    %
    %   Properties:
    %     MP_n(i) is real, ∈ [-1, +1]
    %     MP_n(i) ≠ 0  → nonzero particle-hole coherence at site i
    %     For a true Majorana zero mode: u_n = v_n*, so MP_n(i) = 2|u_n(i)|²
    %     MP_n(i) = 0 everywhere → equal superposition of e and h (not Majorana)
    %
    %   Note: C_n and MP_n are complementary.
    %     C_n measures electron-hole *population* difference.
    %     MP_n measures electron-hole *coherence* (off-diagonal in Nambu space).
    %
    %   ================================================================
    %   NONLOCAL MAJORANA POLARIZATION  MP_nl(i,j)
    %   ================================================================
    %   Two-site correlation between electron at i and hole at j:
    %
    %     MP_nl_n(i,j) = u_n(i)·v_n*(j) + u_n*(i)·v_n(j)
    %                  = 2 Re[ u_n(i)·v_n*(j) ]           (spinless)
    %
    %     MP_nl_n(i,j) = 2 Σ_σ Re[ u_{nσ}(i)·v_{nσ}*(j) ]  (spinful)
    %
    %   Properties:
    %     MP_nl(i,i) = MP_n(i)       (diagonal = local MP)
    %     MP_nl(i,j) = MP_nl(j,i)   (symmetric, real)
    %     For a MZM pair: MP_nl(i_L, i_R) nonzero even if MP(i) ≈ 0 locally
    %
    %   ================================================================
    %   NONLOCAL LEFT-RIGHT POLARIZATION  MP_LR
    %   ================================================================
    %   Product of left and right half-sums:
    %
    %     MP_L  = Σ_{i ≤ N/2} MP_n(i)   (real scalar, left half)
    %     MP_R  = Σ_{i > N/2} MP_n(i)   (real scalar, right half)
    %     MP_LR = MP_L · MP_R            (real scalar, signed product)
    %
    %   MP_LR < 0  : MP_L and MP_R have opposite signs
    %                → topological MZM (left and right modes out of phase)
    %   MP_LR > 0  : same sign on both ends → trivial in-gap state
    %   MP_LR = 0  : no nonlocal correlation
    %
    %   Ref: Awoga O.A., Cayao J., PRB (2023).
    %        Kumar H., Cayao J., Awoga O.A., PRB (2026).
    %
    %   ================================================================
    %   NONLOCAL CORRELATION MATRIX  MP_nl(i,j)
    %   ================================================================
    %   Full N×N matrix of cross-site electron-hole coherence.
    %   Cross-block sum (L↔R):
    %     NMP = Σ_{i≤N/2, j>N/2} |MP_nl(i,j)|
    %   NMP large → strong nonlocal Majorana correlations between ends.
    %
    %   ================================================================
    %   Constructor
    %   ================================================================
    %   Takes eigenvectors directly — no System object required.
    %
    %   mp = MajoranaPolarization(V, N)       % without eigenvalues
    %   mp = MajoranaPolarization(V, N, E)    % with eigenvalues (enables auto tol)
    %
    %   Inputs:
    %     V : Hdim × Nstates  (Hdim = 2N spinless, or 4N spinful BdG)
    %     N : number of physical sites
    %     E : Nstates×1 eigenvalues (optional)
    %
    %   Static convenience:
    %   mp = MajoranaPolarization.fromSystem(sys)
    %
    %   ================================================================
    %   Usage
    %   ================================================================
    %     % Charge C_n(i) = |u|² - |v|²
    %     C  = mp.charge()                % Nstates × N
    %     C  = mp.charge(idx)             % selected states
    %     [C, E, idx] = mp.charge_edge()  % in-gap states only
    %
    %     % Local MP_n(i) = 2 Re[u·v*]
    %     MP = mp.local()                 % Nstates × N
    %     MP = mp.local(idx)
    %     [MP, E, idx] = mp.local_edge()
    %
    %     % Left/right sums and MP_LR = MP_L·conj(MP_R)
    %     [MP_L, MP_R, MP_LR] = mp.left_right(idx)
    %     [P, M] = mp.total(idx)          % P=Σ|MP|, M=P/N
    %
    %     % Nonlocal matrix MP_nl(i,j) = 2 Re[u(i)·v*(j)]
    %     MP_nl = mp.nonlocal(idx)              % N×N, averaged over idx
    %     [MP_nl, E, idx] = mp.nonlocal_edge()  % in-gap states
    %
    %     % Cross-block nonlocal total and MP_LR
    %     [NMP, MP_LR, MP_nl] = mp.nonlocal_total(idx)
    %
    %     % Nonlocal profile: MP_nl(:, j_ref) vs i
    %     [prof, x] = mp.nonlocal_profile(j_ref, idx)
    %
    %     % Electron-hole ratio
    %     r = mp.electron_hole_ratio(idx)
    %
    %     disp(mp)   % full summary: C, MP, MP_LR, NMP for in-gap states
    %
    %  Citation
    %  --------
    %  Local Majorana polarization:
    %    Sticlet D. et al., PRL 108, 096802 (2012). [original definition]
    %    DeGottardi W. et al., PRL 110, 146404 (2013). [spatial structure]
    %  Nonlocal left-right MP_LR = MP_L·conj(MP_R):
    %    Awoga O.A., Cayao J., Phys. Rev. B (2023).
    %    Kumar H., Cayao J., Awoga O.A., Phys. Rev. B (2026).
    %  Nonlocal correlation MP_nl(i,j):
    %    Sedlmayr N. et al., PRB 93, 245109 (2016).
    %    Bena C. et al., PRB 86, 174525 (2012).
    %  Topological Majorana modes:
    %    Kitaev A.Yu., Phys. Uspekhi 44, 131 (2001).
    %  Experimental context:
    %    Mourik V. et al., Science 336, 1003 (2012).
    %    Nadj-Perge S. et al., Science 346, 602 (2014).

    properties (SetAccess = private)
        V       (:,:)          % Hdim × Nstates eigenvector matrix
        E       (:,1) double   % Nstates × 1 eigenvalues (may be empty)
        N       (1,1) double   % number of physical sites
        Hdim    (1,1) double   % Hamiltonian dimension (2N or 4N)
        spinful (1,1) logical  % true = 4N spinful BdG, false = 2N spinless
    end

    methods

        % ============================================================
        % Constructor
        % ============================================================
        function obj = MajoranaPolarization(V, N, E)
            validateattributes(V, {'double'}, {'2d'}, 'MajoranaPolarization', 'V');
            validateattributes(N, {'double'}, {'scalar','positive','integer'}, ...
                'MajoranaPolarization', 'N');
            Hdim = size(V, 1);
            if Hdim == 2*N
                spinful = false;
            elseif Hdim == 4*N
                spinful = true;
            elseif Hdim > 2*N
                % Mixed heterostructure: Hdim is not a clean multiple of N
                % (e.g. SC substrate 4*Ns2d + normal chain Nc).
                % Treat as spinful and use Hdim/2 as the block stride.
                spinful = true;
            else
                error('UQMaT:MP:badSize', ...
                    'V has %d rows but N=%d. Expected 2N=%d (spinless) or 4N=%d (spinful).', ...
                    Hdim, N, 2*N, 4*N);
            end
            obj.V       = V;
            obj.N       = N;
            obj.Hdim    = Hdim;
            obj.spinful = spinful;
            obj.E       = [];
            if nargin >= 3 && ~isempty(E)
                obj.E = E(:);
            end
        end

        % ============================================================
        % CHARGE  C_n(i) = |u_n(i)|² - |v_n(i)|²
        % ============================================================
        function C = charge(obj, state_idx)
            % Local charge per site: C_n(i) = |u_n(i)|² - |v_n(i)|²
            % Returns Nsel × N matrix.
            %   C = mp.charge()          % all states
            %   C = mp.charge(idx)       % selected states
            if nargin < 2 || isempty(state_idx)
                state_idx = 1:size(obj.V, 2);
            end
            Nst = numel(state_idx);
            C   = zeros(Nst, obj.N);
            for k = 1:Nst
                [u2, v2] = obj.split_weights(obj.V(:, state_idx(k)));
                C(k,:)   = u2 - v2;
            end
        end

        function [C, E_edge, idx] = charge_edge(obj, tol)
            % Charge of in-gap states (auto threshold).
            [~, E_edge, idx] = obj.local_edge(nargin>1, tol);
            if isempty(idx), C = []; return; end
            C = obj.charge(idx);
        end

        % ============================================================
        % LOCAL MP  MP_n(i) = 2 Re[Σ_σ u_{nσ}(i)·v_{nσ}*(i)]
        % ============================================================
        function MP = local(obj, state_idx)
            % Local Majorana polarization per site.
            %   Spinless: MP_n(i) = 2 Re[ u_n(i)·v_n*(i) ]
            %   Spinful:  MP_n(i) = 2 Σ_σ Re[ u_{nσ}(i)·v_{nσ}*(i) ]
            % Returns Nsel × N real matrix.
            if nargin < 2 || isempty(state_idx)
                state_idx = 1:size(obj.V, 2);
            end
            Nst = numel(state_idx);
            MP  = zeros(Nst, obj.N);
            for k = 1:Nst
                [u, v] = obj.split_amplitudes(obj.V(:, state_idx(k)));
                % MP(i) = 2 Re[ Σ_σ u_σ(i)·v_σ*(i) ]
                % u and v are N×Nspin matrices; sum over spin then take real
                MP(k,:) = 2*real(sum(u .* conj(v), 2))';
            end
        end

        function [MP, E_edge, idx] = local_edge(obj, tol)
            % Local MP of in-gap states (auto threshold on eigenvalues).
            %   [MP, E_edge, idx] = mp.local_edge()
            %   [MP, E_edge, idx] = mp.local_edge(tol)
            if isempty(obj.E)
                error('UQMaT:MP:noE', ...
                    'Eigenvalues not provided. Pass E to constructor or use local(idx).');
            end
            if nargin < 2 || isempty(tol)
                tol = obj.auto_tol();
            end
            idx    = find(abs(obj.E) < tol);
            E_edge = obj.E(idx);
            if isempty(idx)
                warning('UQMaT:MP:noEdge', ...
                    'No in-gap states found with |E| < %.4g. Try increasing tol.', tol);
                MP = []; return
            end
            MP = obj.local(idx);
        end

        % ============================================================
        % TOTAL  P_n = Σ|MP|,  M_n = P_n/N
        % ============================================================
        function [P, M] = total(obj, state_idx)
            % P(n) = Σ_i |MP_n(i)|    M(n) = P(n)/N ∈ [0,1]
            % M → 0: strong Majorana character
            % M → 1: fermionic (pure electron or hole)
            if nargin < 2 || isempty(state_idx), MP = obj.local();
            else, MP = obj.local(state_idx); end
            P = sum(abs(MP), 2);
            M = P / obj.N;
        end

        % ============================================================
        % LEFT/RIGHT  MP_L, MP_R, MP_LR = MP_L·conj(MP_R)
        % ============================================================
        function [MP_L, MP_R, MP_LR] = left_right(obj, state_idx)
            % Sum local MP over left and right halves.
            %
            %   MP_L(n)  = Σ_{i ≤ N/2} MP_n(i)   (real scalar per state)
            %   MP_R(n)  = Σ_{i > N/2} MP_n(i)   (real scalar per state)
            %   MP_LR(n) = MP_L(n) · MP_R(n)      (real, signed product)
            %
            %   MP_LR < 0 : opposite signs on left and right → topological MZM
            %   MP_LR > 0 : same sign on both ends → trivial
            %   MP_LR = 0 : no nonlocal correlation
            %
            % Ref: Awoga & Cayao PRB (2023); Kumar, Cayao & Awoga PRB (2026).
            if nargin < 2 || isempty(state_idx)
                state_idx = 1:size(obj.V, 2);
            end
            MP    = obj.local(state_idx);
            mid   = floor(obj.N / 2);
            MP_L  = sum(MP(:, 1:mid),     2);   % real
            MP_R  = sum(MP(:, mid+1:end), 2);   % real
            MP_LR = MP_L .* MP_R;               % real, signed
        end

        % ============================================================
        % PROFILE  |MP(i)| averaged over edge states
        % ============================================================
        function [mp_avg, x, MP_all] = profile(obj, tol)
            % Site-resolved average |MP(i)| over in-gap states.
            if nargin < 2, [MP_all] = obj.local_edge();
            else,          [MP_all] = obj.local_edge(tol); end
            if isempty(MP_all), mp_avg = []; x = []; return; end
            mp_avg = mean(abs(MP_all), 1)';
            x      = (1:obj.N)';
        end

        % ============================================================
        % NONLOCAL MATRIX  MP_nl(i,j) = 2 Re[Σ_σ u_{nσ}(i)·v_{nσ}*(j)]
        % ============================================================
        function MP_nl = nonlocal(obj, state_idx)
            % Nonlocal Majorana polarization: N×N real matrix.
            %
            %   Spinless: MP_nl_n(i,j) = 2 Re[ u_n(i)·v_n*(j) ]
            %   Spinful:  MP_nl_n(i,j) = 2 Σ_σ Re[ u_{nσ}(i)·v_{nσ}*(j) ]
            %
            % Properties:
            %   MP_nl(i,i) = MP_n(i)    (diagonal = local MP)
            %   MP_nl(i,j) = MP_nl(j,i) (symmetric, real)
            %
            % Returns N×N matrix averaged over selected states.
            %   MP_nl = mp.nonlocal(n)         % single state
            %   MP_nl = mp.nonlocal([n1,n2])   % averaged over n1,n2
            %   MP_nl = mp.nonlocal()          % all states averaged
            if nargin < 2 || isempty(state_idx)
                state_idx = 1:size(obj.V, 2);
            end
            N     = obj.N;
            MP_nl = zeros(N);
            for k = 1:numel(state_idx)
                [u, v] = obj.split_amplitudes(obj.V(:, state_idx(k)));
                % MP_nl(i,j) = 2 Re[ Σ_σ u_σ(i)·v_σ*(j) ]
                % u, v: N×Nspin;  u*v' is N×N cross-site matrix
                cross = zeros(N);
                for sp = 1:size(u,2)
                    cross = cross + u(:,sp) * conj(v(:,sp))';
                end
                MP_nl = MP_nl + 2*real(cross);
            end
            MP_nl = MP_nl / numel(state_idx);
        end

        function [MP_nl, E_edge, idx] = nonlocal_edge(obj, tol)
            % Nonlocal MP matrix of in-gap states (averaged).
            if nargin < 2, [~, E_edge, idx] = obj.local_edge();
            else,          [~, E_edge, idx] = obj.local_edge(tol); end
            if isempty(idx), MP_nl = []; return; end
            MP_nl = obj.nonlocal(idx);
        end

        % ============================================================
        % NONLOCAL TOTAL  NMP + MP_LR from nonlocal matrix
        % ============================================================
        function [NMP, MP_LR, MP_nl] = nonlocal_total(obj, state_idx)
            % Cross-block nonlocal polarization from the MP_nl(i,j) matrix.
            %
            %   MP_nl : N×N real matrix (averaged over state_idx)
            %   NMP   = Σ_{i≤N/2, j>N/2} MP_nl(i,j)   (signed sum of L-R block)
            %   MP_LR = Sum_LR · Sum_RL                 (real, signed product)
            %         where Sum_LR = Σ_{i≤N/2, j>N/2} MP_nl(i,j)
            %               Sum_RL = Σ_{i>N/2, j≤N/2} MP_nl(i,j)
            %
            %   MP_LR < 0 : L and R nonlocal blocks have opposite sign → topological
            %   MP_LR > 0 : same sign → trivial
            %
            % Ref: Awoga & Cayao PRB (2023); Kumar, Cayao & Awoga PRB (2026).
            if nargin < 2 || isempty(state_idx)
                state_idx = 1:size(obj.V, 2);
            end
            MP_nl    = obj.nonlocal(state_idx);
            N        = obj.N;
            mid      = floor(N/2);
            LR_block = MP_nl(1:mid,     mid+1:end);   % real
            RL_block = MP_nl(mid+1:end, 1:mid);       % real
            Sum_LR   = sum(LR_block(:));               % real scalar
            Sum_RL   = sum(RL_block(:));               % real scalar
            NMP      = Sum_LR;                         % signed sum
            MP_LR    = Sum_LR * Sum_RL;               % real, signed product
        end

        % ============================================================
        % NONLOCAL PROFILE  MP_nl(:, j_ref) vs i
        % ============================================================
        function [nl_prof, x] = nonlocal_profile(obj, j_ref, state_idx)
            % Site-resolved nonlocal MP fixing reference site j_ref.
            % Shows nonlocal electron-hole coherence between j_ref and all i.
            %   [prof, x] = mp.nonlocal_profile(j_ref)
            %   [prof, x] = mp.nonlocal_profile(j_ref, idx)
            %
            % Typical use: j_ref=1 or j_ref=N to probe one end vs all sites.
            if nargin < 2 || isempty(j_ref), j_ref = obj.N; end
            if nargin < 3 || isempty(state_idx)
                [~, ~, state_idx] = obj.local_edge();
                if isempty(state_idx), state_idx = 1:size(obj.V,2); end
            end
            MP_nl   = obj.nonlocal(state_idx);
            nl_prof = MP_nl(:, j_ref);
            x       = (1:obj.N)';
        end

        % ============================================================
        % ELECTRON-HOLE RATIO  |u|²/|v|²
        % ============================================================
        function r = electron_hole_ratio(obj, state_idx)
            % r(n) = Σ_i|u_n(i)|² / Σ_i|v_n(i)|²
            % r = 1 for an ideal MZM (equal electron and hole weight)
            if nargin < 2, state_idx = 1:size(obj.V,2); end
            r = zeros(numel(state_idx), 1);
            for k = 1:numel(state_idx)
                [u2, v2] = obj.split_weights(obj.V(:, state_idx(k)));
                tot_u = sum(u2); tot_v = sum(v2);
                if tot_v < 1e-15, r(k) = Inf;
                else,              r(k) = tot_u/tot_v; end
            end
        end

        % ============================================================
        % disp
        % ============================================================
        function disp(obj)
            fprintf('\n  UQMaT MajoranaPolarization\n');
            fprintf('  N=%d  Hdim=%d  %s BdG\n', ...
                obj.N, obj.Hdim, iif(obj.spinful,'spinful','spinless'));
            if isempty(obj.E)
                fprintf('  No eigenvalues stored — use mp.local(idx), mp.charge(idx).\n\n');
                return
            end
            try
                [~, E_edge, idx] = obj.local_edge();
            catch
                fprintf('  Cannot determine in-gap states.\n\n'); return
            end
            if isempty(idx)
                fprintf('  No in-gap states found.\n\n'); return
            end
            fprintf('  In-gap states: %d  E = [', numel(idx));
            fprintf(' %.4g', E_edge); fprintf(' ]\n\n');

            fprintf('  %-6s  %-8s  %-8s  %-8s  %-8s  %-10s\n', ...
                'State', 'P', 'M', 'MP_L', 'MP_R', 'MP_LR');
            fprintf('  %s\n', repmat('-',1,58));
            [P, M]              = obj.total(idx);
            [MP_L, MP_R, MP_LR] = obj.left_right(idx);
            for k = 1:numel(idx)
                fprintf('  %-6d  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-10.4f\n', ...
                    idx(k), P(k), M(k), MP_L(k), MP_R(k), MP_LR(k));
            end
            fprintf('\n  M→0: Majorana-like    M→1: fermionic\n');
            fprintf('  MP_LR<0: opposite L/R polarization → topological MZM\n');
            fprintf('  MP_LR>0: same-sign polarization → trivial\n');

            % Nonlocal summary
            [NMP, MP_LR_nl] = obj.nonlocal_total(idx);
            fprintf('\n  Nonlocal (L↔R block):\n');
            fprintf('  NMP = %.4f   MP_LR_nl = %.4f\n', NMP, MP_LR_nl);
            fprintf('\n');
        end

    end   % public methods

    methods (Access = private)

        function [u2, v2] = split_weights(obj, psi)
            % |u(i)|² and |v(i)|² per site (summed over spin).
            % For mixed heterostructures Hdim may not equal 2N or 4N exactly;
            % use Hdim/2 as the electron/hole block boundary.
            N    = obj.N;
            Hdim = obj.Hdim;
            half = floor(Hdim / 2);   % electron block: 1:half, hole block: half+1:end
            if obj.spinful
                % Try standard 4N layout first
                if Hdim == 4*N
                    u2 = abs(psi(1:N)).^2     + abs(psi(N+1:2*N)).^2;
                    v2 = abs(psi(2*N+1:3*N)).^2 + abs(psi(3*N+1:4*N)).^2;
                else
                    % Mixed: split at Hdim/2, sum last N rows of each half
                    e_block = psi(1:half);
                    h_block = psi(half+1:end);
                    u2 = abs(e_block(end-N+1:end)).^2;
                    v2 = abs(h_block(end-N+1:end)).^2;
                end
            else
                u2 = abs(psi(1:N)).^2;
                v2 = abs(psi(N+1:2*N)).^2;
            end
        end

        function [u, v] = split_amplitudes(obj, psi)
            % Complex amplitudes for MP and MP_nl computation.
            N    = obj.N;
            Hdim = obj.Hdim;
            half = floor(Hdim / 2);
            if obj.spinful
                if Hdim == 4*N
                    u = [psi(1:N), psi(N+1:2*N)];         % N×2
                    v = [psi(2*N+1:3*N), psi(3*N+1:4*N)]; % N×2
                else
                    % Mixed: take last N rows of each half as representative
                    e_block = psi(1:half);
                    h_block = psi(half+1:end);
                    u = e_block(end-N+1:end);   % N×1
                    v = h_block(end-N+1:end);   % N×1
                end
            else
                u = psi(1:N);
                v = psi(N+1:2*N);
            end
        end

        function tol = auto_tol(obj)
            pos = obj.E(obj.E >  1e-10);
            neg = obj.E(obj.E < -1e-10);
            if ~isempty(pos) && ~isempty(neg)
                tol = min(pos(1), abs(neg(end))) * 0.5;
            else
                tol = 0.05 * (obj.E(end) - obj.E(1));
            end
        end

    end   % private methods

    methods (Static)

        function mp = fromSystem(sys)
            % Convenience constructor from a UQMaT.System.
            %   mp = MajoranaPolarization.fromSystem(sys)
            if ~isa(sys, 'UQMaT.System')
                error('UQMaT:MP:badInput', 'Argument must be a UQMaT.System.');
            end
            if ~sys.params.sc && ~ismember(sys.params.system,{'kitaev','hnkitaev'})
                error('UQMaT:MP:notBdG', ...
                    'MajoranaPolarization requires sc=true or system=kitaev/hnkitaev. Got "%s".', ...
                    sys.params.system);
            end
            mp = UQMaT.Observables.Topological.RealSpace.MajoranaPolarization(...
                sys.V, sys.Ndim, sys.E);
        end

    end   % static methods
end

function s = iif(cond, a, b)
    if cond, s = a; else, s = b; end
end
