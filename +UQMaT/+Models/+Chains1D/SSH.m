classdef SSH < UQMaT.Models.Base
    % UQMaT.Models.Chains1D.SSH
    %   Su-Schrieffer-Heeger chain (dimerised 1D lattice).
    %   2 sites per unit cell: A and B, alternating hoppings v (intra) and w (inter).
    %
    %   Physical sites: Ns = 2*N  (N unit cells)
    %   Basis: [A_1, B_1, A_2, B_2, ..., A_N, B_N]
    %
    %   H sizes — determined automatically from parameters:
    %     sc=false, spin=false  →  2N × 2N     normal sublattice
    %     sc=false, spin=true   →  4N × 4N     spinful normal
    %     sc=true               →  8N × 8N     BdG Ψ=(c↑,c↓,c↑†,c↓†)ᵀ per site
    %
    %   Topological phase (normal): w > v  (Zak phase = π, edge states)
    %   Trivial phase (normal):     w < v  (Zak phase = 0)
    %
    %   With sc=true: proximity-induced SC on the SSH chain.
    %   Error raised if sc=true but all pairing amplitudes are zero.

    %  Citation
    %  --------
    %  Su W.P., Schrieffer J.R., Heeger A.J., PRL 42, 1698 (1979).
    %  Su W.P., Schrieffer J.R., Heeger A.J., PRB 22, 2099 (1980).
    %  Asboth J.K., Oroszlany L., Palyi A., "A Short Course on Topological
    %    Insulators", Springer (2016). [SSH model, Zak phase, edge states]
    methods
        function obj = SSH(p)
            obj.params = p;
            obj.Ndim   = 2 * p.N;   % Ns = 2N physical sites
            obj.H      = obj.buildH();
        end

        function s = description(obj)
            p  = obj.params;
            Ns = 2 * p.N;
            if p.w > p.v,     ph = 'topological';
            elseif p.w < p.v, ph = 'trivial';
            else,              ph = 'critical';
            end
            if p.sc
                s = sprintf('SSH BdG 8×%d (%d cells) bc=%s v=%.3g w=%.3g [%s] Δ_s=%.3g', ...
                    p.N, p.N, upper(p.bc), p.v, p.w, ph, mean(p.Delta_s_vec));
            elseif p.spin
                s = sprintf('SSH spinful 4×%d (%d cells) bc=%s v=%.3g w=%.3g [%s]', ...
                    p.N, p.N, upper(p.bc), p.v, p.w, ph);
            else
                s = sprintf('SSH %d cells (%d sites) bc=%s v=%.3g w=%.3g [%s]', ...
                    p.N, Ns, upper(p.bc), p.v, p.w, ph);
            end
        end

        function H = buildH(obj)
            p   = obj.params;
            N   = p.N;
            Ns  = 2 * N;   % physical sites
            bdg = UQMaT.SC.BdG;

            % Validate: sc=true requires at least one nonzero amplitude
            if p.sc && all(p.Delta_s_vec == 0) && p.Delta_t_val == 0
                error('UQMaT:SSH:noAmplitude', ...
                    ['sc=true but all pairing amplitudes are zero.\n' ...
                     'Set Delta_s, dvec, or Delta_t to a nonzero value.\n' ...
                     'Example: UQMaT.Params(''system'',''ssh'',''Delta_s'',0.3)']);
            end

            % Normal sublattice Hamiltonian (Ns × Ns)
            H0 = obj.buildNormal(N, p);

            if p.sc
                % ---- BdG: 4Ns × 4Ns ----
                H_N = kron(eye(2), H0);
                if any(p.B ~= 0)
                    H_N = H_N + bdg.zeeman(p.B, Ns);
                end
                if p.lambda_R ~= 0
                    Tr  = obj.rashbaHop(N, p);
                    H_N = H_N + bdg.rashba(p.lambda_R, Tr);
                end
                D = obj.pairingBlock(Ns, p);
                H = bdg.assemble(H_N, D);

            elseif p.spin
                % ---- Spinful normal: 2Ns × 2Ns ----
                H = kron(eye(2), H0);
                if any(p.B ~= 0)
                    H = H + bdg.zeeman(p.B, Ns);
                end
                if p.lambda_R ~= 0
                    Tr = obj.rashbaHop(N, p);
                    H  = H + bdg.rashba(p.lambda_R, Tr);
                end

            else
                % ---- Normal sublattice: Ns × Ns ----
                H = H0;
            end
        end
    end

    methods (Access = private)

        % ---- Normal SSH Hamiltonian (Ns × Ns) ----------------------
        function H0 = buildNormal(obj, N, p)
            dim = 2 * N;
            H0  = zeros(dim);
            for n = 1:N
                iA = 2*n-1;  iB = 2*n;
                H0(iA,iA) = p.mu + p.onsite_val;
                H0(iB,iB) = p.mu + p.onsite_val;
                H0(iA,iB) = -p.v;  H0(iB,iA) = -p.v;   % intra-cell
                if n < N
                    iA2 = 2*n+1;
                    H0(iB,iA2) = -p.w;  H0(iA2,iB) = -p.w;  % inter-cell
                end
            end
            if strcmp(p.bc,'pbc')
                H0(2*N,1) = -p.w;
                H0(1,2*N) = -p.w;
            end
            % On-site disorder via Disorder class
            H0 = H0 + diag(UQMaT.Disorder.Disorder.onsite(p, dim));
        end

        % ---- Rashba antisymmetric matrix on SSH bonds (Ns × Ns) ----
        function A = rashbaHop(~, N, p)
            % Antisymmetric nearest-neighbour matrix on the SSH chain.
            % A(i,i+1) = +1, A(i+1,i) = -1 for each bond.
            Ns = 2*N;
            A  = diag(ones(Ns-1,1),1) - diag(ones(Ns-1,1),-1);
            if strcmp(p.bc,'pbc')
                A(Ns,1) = -1;   A(1,Ns) = +1;
            end
        end

        % ---- Pairing block (2Ns × 2Ns) ----------------------------
        function D = pairingBlock(~, Ns, p)
            bdg    = UQMaT.SC.BdG;
            % Delta_s_vec was expanded to length p.N in Params.
            % Here Ns = 2*p.N, so expand by repeating (uniform per cell).
            Ds = p.Delta_s_vec;
            if numel(Ds) < Ns
                Ds = repmat(Ds(1), Ns, 1);
            end
            switch p.pairing_mode
                case 'singlet'
                    D = bdg.singlet(Ds(1:Ns), Ns);
                case 'triplet'
                    D = bdg.triplet(zeros(1,3), Ns, p.Delta_t_val, p.bc);
                case 'mixed'
                    D = bdg.singlet(Ds(1:Ns), Ns) + ...
                        bdg.triplet(zeros(1,3), Ns, p.Delta_t_val, p.bc);
                otherwise
                    error('UQMaT:SSH:badPairing', ...
                        'Unknown pairing_mode: "%s"', p.pairing_mode);
            end
        end

    end
end
