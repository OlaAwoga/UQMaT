classdef Hamiltonian < UQMaT.Models.Base
    % UQMaT.Models.Lattice.Nambu.Hamiltonian
    %   Full 4-component Nambu BdG Hamiltonian for a 1D spinful SC chain.
    %
    %   Nambu spinor:  Ψ = (c↑, c↓, c↑†, c↓†)ᵀ
    %
    %   Block structure (4N × 4N):
    %
    %     H_BdG = [  H_N(k)       Δ̂        ]
    %             [ -Δ̂†        -H_N(-k)*  ]
    %
    %     H_N = h_k ⊗ I₂  +  B⃗·σ⃗  +  H_Rashba     (2N × 2N)
    %
    %     Δ̂   = Δ_s·(iσ_y)⊗I_N                     (singlet, s-wave)
    %           + (d⃗·σ⃗)(iσ_y)⊗T_N                  (triplet, p-wave)
    %
    %   Supported pairing channels (set via Params.pairing):
    %     'singlet'   — s-wave only            Δ̂ = Δ_s·iσ_y
    %     'triplet'   — p-wave only             Δ̂ = (d⃗·σ⃗)(iσ_y)·T
    %     'mixed'     — singlet + triplet       Δ̂ = both
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params( ...
    %               'N', 80, 'mu', -1, 't', 1, 'bc', 'obc', ...
    %               'pairing', 'singlet', 'Delta_s', 0.3, ...
    %               'Bz', 0.2, 'lambda_R', 0.4);
    %     mdl = UQMaT.Models.Lattice.Nambu.Hamiltonian(p);
    %     sys = UQMaT.System(mdl);
    %
    %   For triplet / mixed:
    %     p   = UQMaT.Params('pairing','triplet','dvec',[0 0 0.5], ...
    %               'N',80,'mu',-1,'t',1);
    %     p   = UQMaT.Params('pairing','mixed', ...
    %               'Delta_s',0.2,'dvec',[0.3 0 0],'N',80);
    %
    %   Index convention (see UQMaT.Models.Lattice.Nambu.Basis for full docs):
    %     rows  1 : N      →  c_{i↑}     (electron spin-up)
    %     rows  N+1 : 2N   →  c_{i↓}     (electron spin-down)
    %     rows  2N+1 : 3N  →  c†_{i↑}   (hole spin-up)
    %     rows  3N+1 : 4N  →  c†_{i↓}   (hole spin-down)
    %
    %   Symmetry checks (call after construction):
    %     UQMaT.Models.Lattice.Nambu.Basis.ph_sym(mdl.H)
    %     UQMaT.Models.Lattice.Nambu.Basis.tr_sym(mdl.H)

    methods
        function obj = Hamiltonian(p)
            if ~isa(p, 'UQMaT.Params')
                error('UQMaT:Nambu:Hamiltonian:badInput', ...
                    'Argument must be a UQMaT.Params object.');
            end
            obj.params = p;
            obj.Ndim   = obj.ndim();
            obj.H      = obj.buildH();
        end

        % ---- Base interface --------------------------------------------
        function n = ndim(obj)
            n = obj.params.N;   % physical sites; H is 4N×4N
        end

        function s = family(~)
            s = 'lattice';
        end

        function s = description(obj)
            p = obj.params;
            s = sprintf( ...
                'Nambu BdG 4×%d | pairing=%s | Δ_s=%.3g | d=[%.3g,%.3g,%.3g] | Bz=%.3g λ_R=%.3g', ...
                p.N, p.pairing_mode, p.sc_Delta_s_s, 0, 0, 0, ...
                p.B(3), p.lambda_R);
        end

        % ---- Core build ------------------------------------------------
        function H = buildH(obj)
            p  = obj.params;
            N  = p.N;
            nb = UQMaT.Models.Lattice.Nambu.Basis;

            % --- 1. Scalar hopping + disorder (N×N) ---
            h_k = obj.hopMatrix(N, p.t, p.t2, p.bc) - p.mu * eye(N);
            h_k = h_k + obj.disorderDiag(N);

            % --- 2. Normal block H_N (2N×2N) ---
            %   H_N = h_k ⊗ I₂  +  Zeeman  +  Rashba
            H_N = kron(eye(2), h_k);                          % spin-diagonal hopping

            if p.B(1) ~= 0 || p.B(2) ~= 0 || p.B(3) ~= 0
                H_N = H_N + nb.zeeman(p.B(1), p.B(2), p.B(3), N);
            end

            if p.lambda_R ~= 0
                T_hop = obj.hopMatrix(N, p.lambda_R, 0, p.bc);
                H_N   = H_N + nb.rashba(p.lambda_R, T_hop);
            end

            % --- 3. Pairing block Δ̂ (2N×2N) ---
            Delta_hat = obj.buildPairing(N, p);

            % --- 4. Assemble 4N×4N BdG ---
            H = nb.assemble(H_N, Delta_hat);
        end

        % ---- Display ---------------------------------------------------
        function disp(obj)
            p = obj.params;
            N = p.N;
            fprintf('\n  UQMaT.Models.Lattice.Nambu.Hamiltonian\n');
            fprintf('  Basis  : Ψ = (c↑, c↓, c↑†, c↓†)ᵀ\n');
            fprintf('  Size   : %d × %d  (4N, N=%d)\n', 4*N, 4*N, N);
            fprintf('  BC     : %s\n', upper(p.bc));
            fprintf('  h_k    : t=%.4g, μ=%.4g, t2=%.4g\n', p.t, p.mu, p.t2);
            fprintf('  Zeeman : Bx=%.4g  By=%.4g  Bz=%.4g\n', p.B(1), p.B(2), p.B(3));
            fprintf('  Rashba : λ_R=%.4g\n', p.lambda_R);
            fprintf('  Pairing: %s\n', p.pairing_mode);
            switch p.pairing_mode
                case {'singlet','mixed'}
                    fprintf('    Δ_s = %.4g (s-wave singlet)\n', p.sc_Delta_s_s);
            end
            switch p.pairing_mode
                case {'triplet','mixed'}
                    fprintf('    d⃗  = [%.4g, %.4g, %.4g] (triplet)\n', ...
                        0, 0, 0);
            end
            if ~strcmp(p.disorder,'none')
                fprintf('  Disorder: %s, W=%.4g\n', p.disorder, p.W);
            end
            fprintf('\n');
            % Auto symmetry check
            UQMaT.Models.Lattice.Nambu.Basis.ph_sym(obj.H);
            UQMaT.Models.Lattice.Nambu.Basis.tr_sym(obj.H);
        end
    end

    % ====================================================================
    methods (Access = private)

        % ---- Pairing block Δ̂ (2N × 2N) --------------------------------
        function D = buildPairing(obj, N, p)
            nb = UQMaT.Models.Lattice.Nambu.Basis;
            D  = zeros(2*N);

            switch p.pairing_mode
                case 'none'
                    % Normal state — no pairing

                case 'singlet'
                    % Δ̂ = Δ_s · iσ_y ⊗ diag(Delta_s_vec)
                    D = nb.singlet_pairing(p.sc_Delta_s_s_vec, N);

                case 'triplet'
                    % Δ̂ = (d⃗·σ⃗)(iσ_y) ⊗ T_hop
                    T_hop = obj.hopMatrix(N, 1.0, 0, p.bc);  % unit bond matrix
                    D     = nb.triplet_pairing([0,0,0], T_hop * diag(p.sc_Delta_s_t_vec), N);

                case 'mixed'
                    % Both singlet and triplet
                    D_s   = nb.singlet_pairing(p.sc_Delta_s_s_vec, N);
                    T_hop = obj.hopMatrix(N, 1.0, 0, p.bc);
                    D_t   = nb.triplet_pairing([0,0,0], T_hop * diag(p.sc_Delta_s_t_vec), N);
                    D     = D_s + D_t;

                otherwise
                    error('UQMaT:Nambu:Hamiltonian:badPairing', ...
                        'pairing must be ''none'',''singlet'',''triplet'', or ''mixed''.');
            end
        end

        % ---- On-site disorder diagonal (N×N) ----------------------------
        function Vd = disorderDiag(obj, N)
            p = obj.params;
            if p.seed > 0, rng(p.seed); end
            switch p.disorder
                case 'none'
                    Vd = zeros(N);
                case 'anderson'
                    Vd = diag(p.W * (rand(N,1) - 0.5));
                case 'binary'
                    Vd = diag(p.W * (2*round(rand(N,1)) - 1) / 2);
                case 'impurity'
                    Vd = zeros(N);
                    for k = 1:numel(p.imp_sites)
                        i = p.imp_sites(k);
                        if i >= 1 && i <= N
                            Vd(i,i) = p.imp_potential(k);
                        end
                    end
                otherwise
                    Vd = zeros(N);
            end
        end

        % ---- 1D hopping matrix (N×N) ------------------------------------
        function T = hopMatrix(~, N, t, t2, bc)
            e1 = -t * ones(N-1, 1);
            T  = diag(e1, 1) + diag(e1, -1);
            if t2 ~= 0
                e2 = -t2 * ones(N-2, 1);
                T  = T + diag(e2, 2) + diag(e2, -2);
            end
            if strcmp(bc, 'pbc')
                T(1,N) = -t;   T(N,1) = -t;
                if t2 ~= 0 && N > 3
                    T(1,N-1)=-t2; T(2,N)=-t2;
                    T(N-1,1)=-t2; T(N,2)=-t2;
                end
            end
        end

    end
end
