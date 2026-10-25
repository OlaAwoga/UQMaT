classdef BdG
    % UQMaT.SC.BdG
    %   BdG Hamiltonian utilities for superconducting systems.
    %
    %   Spinor basis:  Ψ_i = (c_{i↑}, c_{i↓}, c†_{i↑}, c†_{i↓})ᵀ
    %
    %   Row blocks (4N × 4N):
    %     1   : N     →  c_{i↑}    electron spin-up
    %     N+1 : 2N    →  c_{i↓}    electron spin-down
    %     2N+1: 3N    →  c†_{i↑}   hole spin-up
    %     3N+1: 4N    →  c†_{i↓}   hole spin-down
    %
    %   BdG structure:
    %     H_BdG = [  H_N,     D  ]
    %             [ -D†,   -H_N* ]
    %
    %     H_N : 2N×2N normal block  (hopping + Zeeman + Rashba)
    %     D   : 2N×2N pairing block (singlet / triplet / mixed)
    %
    %   Pairing:
    %     Singlet s-wave : D = Delta_s · (iσ_y) ⊗ I_N
    %     Triplet p-wave : D = (d·σ)(iσ_y) ⊗ A_N
    %       where A_N is antisymmetric: A_N(i,i+1)=+Δ_t, A_N(i+1,i)=-Δ_t
    %
    %   Particle-hole symmetry: Ξ H_BdG Ξ⁻¹ = -H_BdG
    %     Ξ = τ_x ⊗ σ_0 · K  (K = complex conjugation)
    %
    %   Usage
    %   -----
    %     sig = UQMaT.SC.BdG.sigma()         % Pauli matrices struct
    %     D   = UQMaT.SC.BdG.singlet(Ds, N)  % singlet pairing block
    %     D   = UQMaT.SC.BdG.triplet(dvec, N, Dt, bc)  % triplet
    %     H   = UQMaT.SC.BdG.assemble(H_N, D)          % 4N×4N BdG
    %     UQMaT.SC.BdG.check_ph(H)           % particle-hole check
    %     UQMaT.SC.BdG.check_tr(H)           % time-reversal check
    %     idx = UQMaT.SC.BdG.elec_up(N)      % row indices electron ↑
    %     idx = UQMaT.SC.BdG.elec_dn(N)      % row indices electron ↓
    %     idx = UQMaT.SC.BdG.hole_up(N)      % row indices hole ↑
    %     idx = UQMaT.SC.BdG.hole_dn(N)      % row indices hole ↓

    %  Citation
    %  --------
    %  Bogoliubov N.N., Nuovo Cimento 7, 794 (1958).
    %  de Gennes P.G., "Superconductivity of Metals and Alloys",
    %    Addison-Wesley (1966).
    %  Zhu J.-X., "Bogoliubov-de Gennes Method and Its Applications",
    %    Springer (2016). [Chapter 1: BdG formalism, basis, and assembly]
    methods (Static)

        % ---- Pauli matrices ----------------------------------------
        function s = sigma()
            s.s0 = eye(2);
            s.sx = [0  1; 1  0];
            s.sy = [0 -1i; 1i 0];
            s.sz = [1  0; 0 -1];
        end

        % ---- Antisymmetric Rashba hopping matrix -------------------
        function A = rashba_hop(N, bc)
            % Build N×N antisymmetric nearest-neighbour matrix:
            %   A(i, i+1) = +1   (forward hop, gains phase)
            %   A(i+1, i) = -1   (backward hop, loses phase)
            % Used as A_hop in rashba(lambda_R, A_hop).
            if nargin < 2, bc = 'obc'; end
            A = diag(ones(N-1,1), 1) - diag(ones(N-1,1), -1);
            if strcmp(bc,'pbc')
                A(N,1) = -1;   % wrap-around: backward
                A(1,N) = +1;   % wrap-around: forward
            end
        end

        % ---- Index helpers -----------------------------------------
        function idx = elec_up(N), idx = 1:N;          end
        function idx = elec_dn(N), idx = (1:N)+N;      end
        function idx = hole_up(N), idx = (1:N)+2*N;    end
        function idx = hole_dn(N), idx = (1:N)+3*N;    end
        function idx = elec_block(N), idx = 1:2*N;     end
        function idx = hole_block(N), idx = (1:2*N)+2*N; end

        % ---- Zeeman term in 2N×2N normal block ---------------------
        function HZ = zeeman(Bx_or_vec, By_or_N, Bz_or_empty, N_or_empty)
            % Accepts two calling conventions:
            %   zeeman(B, N)          B = [Bx, By, Bz],  N = number of sites
            %   zeeman(Bx, By, Bz, N) three scalars + N
            s = UQMaT.SC.BdG.sigma();
            if nargin == 2
                % zeeman(B, N)
                B  = Bx_or_vec(:);
                N  = By_or_N;
                HZ = kron(B(1)*s.sx + B(2)*s.sy + B(3)*s.sz, eye(N));
            else
                % zeeman(Bx, By, Bz, N)
                Bx = Bx_or_vec;  By = By_or_N;  Bz = Bz_or_empty;
                N  = N_or_empty;
                HZ = kron(Bx*s.sx + By*s.sy + Bz*s.sz, eye(N));
            end
        end

        % ---- Rashba SOC in 2N×2N normal block ----------------------
        function HR = rashba(lambda_R, A_hop)
            % HR = i*lambda_R * (σ_y ⊗ A)
            %
            % A_hop MUST be antisymmetric: A(i,i+1) = +1, A(i+1,i) = -1
            % (the finite-difference representation of d/dx)
            %
            % This gives a Hermitian Rashba term:
            %   HR† = conj(i*lambda_R) * kron(sy†, A^T)
            %       = (-i*lambda_R) * kron(sy, -A)
            %       = i*lambda_R * kron(sy, A) = HR ✓
            %
            % NOTE: A_hop must be antisymmetric, not the symmetric hopping T.
            %   Wrong: T(i,i+1) = T(i+1,i) = -t  (symmetric → HR not Hermitian)
            %   Right: A(i,i+1) = +1, A(i+1,i) = -1  (antisymmetric → HR Hermitian)
            s  = UQMaT.SC.BdG.sigma();
            HR = 1i * lambda_R * kron(s.sy, A_hop);
        end

        % ---- Singlet s-wave pairing block (2N×2N) ------------------
        function D = singlet(Delta_vec, N)
            % D = iσ_y ⊗ diag(Delta_vec)
            % iσ_y = [0,1;-1,0] → antisymmetric in spin = singlet
            s = UQMaT.SC.BdG.sigma();
            D = kron(1i*s.sy, diag(Delta_vec(:)));
        end

        % ---- Triplet p-wave pairing block (2N×2N) ------------------
        function D = triplet(dvec, N, Delta_t, bc)
            % D = (d·σ)(iσ_y) ⊗ A_N
            % A_N is antisymmetric: A(i,i+1)=+Dt, A(i+1,i)=-Dt
            % This satisfies fermionic antisymmetry: D(i,j)=-D(j,i)
            if nargin < 4, bc = 'obc'; end
            s  = UQMaT.SC.BdG.sigma();
            dx = dvec(1); dy = dvec(2); dz = dvec(3);
            % (d·σ)(iσ_y)
            dsisy = [dz,        dx-1i*dy;
                     dx+1i*dy, -dz      ];
            % Antisymmetric bond matrix
            A = zeros(N);
            for i = 1:N-1
                A(i,  i+1) = +Delta_t;
                A(i+1,i  ) = -Delta_t;
            end
            if strcmp(bc,'pbc')
                A(N,1) = +Delta_t;
                A(1,N) = -Delta_t;
            end
            D = kron(dsisy, A);
        end

        % ---- Assemble 4N×4N BdG ------------------------------------
        function H = assemble(H_N, D)
            % H_BdG = [  H_N,   D  ]
            %         [  D†,  -H_N* ]
            %
            % Requirements for a valid BdG Hamiltonian:
            %   1. H_N must be Hermitian (H_N† = H_N)
            %   2. D must be antisymmetric in spin×site (D = -D^T, not just D†)
            %
            % With these satisfied:
            %   H_BdG† = [H_N†, (D†)†] = [H_N, D] = upper blocks ✓
            %            [D†, (-H_N*)†]   [D†, -H_N^T] ← Hermitian iff H_N real
            %
            % The correct assembly is:
            %   Lower-left  = +D†        (NOT -D†)
            %   Lower-right = -conj(H_N) (conjugate, NOT transpose)
            %
            % This gives a Hermitian H_BdG when H_N is Hermitian
            % and is PH symmetric: τ_x conj(H_BdG) τ_x = -H_BdG
            H = [ H_N,       D  ;
                  D',   -conj(H_N) ];
        end

        % ---- Particle-hole symmetry check --------------------------
        function [ok, res] = check_ph(H)
            dim = size(H,1);  N2 = dim/2;
            Tx  = [zeros(N2), eye(N2); eye(N2), zeros(N2)];
            res = norm(Tx*conj(H)*Tx + H,'fro') / norm(H,'fro');
            ok  = res < 1e-10;
            if ok
                fprintf('  PH symmetry: OK  (residual=%.2e)\n', res);
            else
                fprintf('  PH symmetry: VIOLATED  (residual=%.2e)\n', res);
            end
        end

        % ---- Time-reversal symmetry check --------------------------
        function [ok, res] = check_tr(H)
            dim = size(H,1);  N4 = dim/4;
            s   = UQMaT.SC.BdG.sigma();
            isy = kron(1i*s.sy, eye(2*N4));
            res = norm(isy*conj(H)*isy' - H,'fro') / (norm(H,'fro')+eps);
            ok  = res < 1e-10;
            if ok
                fprintf('  TR symmetry: OK  (residual=%.2e)\n', res);
            else
                fprintf('  TR symmetry: broken  (residual=%.2e)\n', res);
            end
        end

        % ---- Print basis info --------------------------------------
        function info()
            fprintf('\n  UQMaT SC basis: Psi_i = (c_up, c_dn, c_up+, c_dn+)\n');
            fprintf('  Rows  1   :  N   -> c_{i,up}   electron spin-up\n');
            fprintf('  Rows N+1  : 2N   -> c_{i,dn}   electron spin-down\n');
            fprintf('  Rows 2N+1 : 3N   -> c+_{i,up}  hole spin-up\n');
            fprintf('  Rows 3N+1 : 4N   -> c+_{i,dn}  hole spin-down\n\n');
            fprintf('  H_BdG = [  H_N,   D  ]\n');
            fprintf('          [ -D+, -H_N* ]\n\n');
            fprintf('  Singlet: D = Delta_s*(iσ_y) ⊗ I_N  (k-independent)\n');
            fprintf('  Triplet: D = (d·σ)(iσ_y) ⊗ A_N    (A antisymmetric)\n\n');
        end

    end
end
