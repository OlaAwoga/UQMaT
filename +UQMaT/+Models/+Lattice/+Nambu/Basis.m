classdef Basis
    % UQMaT.Models.Lattice.Nambu.Basis
    %   Nambu basis conventions, Pauli matrices, and BdG assembly utilities.
    %   Shared by all superconducting models in UQMaT.
    %
    % ===================================================================
    %  NAMBU SPINOR
    % ===================================================================
    %
    %   Ψ_i = ( c_{i↑},  c_{i↓},  c†_{i↑},  c†_{i↓} )ᵀ
    %
    %   Row blocks in the 4N×4N (or 4Ns×4Ns) BdG matrix:
    %
    %     rows  1    : N     →  c_{i↑}     electron spin-up
    %     rows  N+1  : 2N    →  c_{i↓}     electron spin-down
    %     rows  2N+1 : 3N    →  c†_{i↑}   hole spin-up
    %     rows  3N+1 : 4N    →  c†_{i↓}   hole spin-down
    %
    % ===================================================================
    %  BdG STRUCTURE
    % ===================================================================
    %
    %   H_BdG = [  H_N,     Δ̂   ]
    %            [ -Δ̂†,  -H_N* ]
    %
    %   H_N (2N×2N) normal block:
    %     H_N = h_k ⊗ I₂  +  B⃗·σ⃗  +  H_Rashba
    %
    %   Δ̂ (2N×2N) pairing block:
    %     Singlet s-wave : Δ̂ = Δ_s · iσ_y ⊗ I_N
    %     Triplet p-wave : Δ̂ = (d⃗·σ⃗)(iσ_y) ⊗ A_N
    %       where A_N is the ANTISYMMETRIC bond hopping matrix
    %       A_N(i,i+1) = +Δ_t,  A_N(i+1,i) = -Δ_t  (p-wave: odd under i↔j)
    %
    % ===================================================================
    %  PARTICLE-HOLE SYMMETRY
    % ===================================================================
    %
    %   Ξ = τ_x ⊗ σ_0 · K    (K = complex conjugation)
    %   Consequence: if (E, Ψ) eigenstate → (-E, Ξ·Ψ*) eigenstate
    %   Spectrum symmetric about E = 0.
    %
    % ===================================================================
    %  TIME-REVERSAL SYMMETRY
    % ===================================================================
    %
    %   Θ = iσ_y ⊗ τ_0 · K
    %   Preserved when B = 0.

    methods (Static)

        % ----------------------------------------------------------------
        % PAULI MATRICES
        % ----------------------------------------------------------------
        function s = sigma()
            % Spin Pauli matrices
            s.s0 = eye(2);
            s.sx = [0  1;  1  0];
            s.sy = [0 -1i; 1i  0];
            s.sz = [1  0;  0 -1];
        end

        function t = tau()
            % Particle-hole Pauli matrices
            t.t0 = eye(2);
            t.tx = [0 1; 1 0];
            t.ty = [0 -1i; 1i 0];
            t.tz = [1 0; 0 -1];
        end

        % ----------------------------------------------------------------
        % NORMAL BLOCK TERMS
        % ----------------------------------------------------------------
        function HZ = zeeman(Bx, By, Bz, N)
            % Zeeman term in 2N×2N normal block
            sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
            HZ  = kron(Bx*sig.sx + By*sig.sy + Bz*sig.sz, eye(N));
        end

        function HR = rashba(lambda_R, T_hop)
            % Rashba SOC in 2N×2N normal block
            % H_R = i·λ_R · σ_y ⊗ T_hop  (nearest-neighbour Rashba)
            sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
            HR  = 1i * lambda_R * kron(sig.sy, T_hop);
        end

        % ----------------------------------------------------------------
        % PAIRING BLOCKS (2N×2N)
        % ----------------------------------------------------------------
        function D = singlet_pairing(Delta_vec, N)
            % Singlet s-wave pairing: Δ̂ = iσ_y ⊗ diag(Delta_vec)
            %
            % iσ_y = [ 0,  1 ]   →   Δ̂ = [ 0,       diag(Δ) ]
            %        [-1,  0 ]            [-diag(Δ),  0       ]
            %
            % This is ANTISYMMETRIC in spin (singlet) and SYMMETRIC on-site.
            sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
            D   = kron(1i * sig.sy, diag(Delta_vec));
        end

        function D = triplet_pairing(dvec, N, Delta_t, bc)
            % Triplet p-wave pairing: Δ̂ = (d⃗·σ⃗)(iσ_y) ⊗ A_N
            %
            % A_N is the ANTISYMMETRIC bond matrix (p-wave requires odd parity):
            %   A_N(i,i+1) = +Delta_t   (forward bond)
            %   A_N(i+1,i) = -Delta_t   (backward bond, antisymmetric)
            %
            % (d⃗·σ⃗)(iσ_y) = [ dz,       dx - i·dy ]
            %                 [ dx+i·dy,   -dz       ]
            %
            % Bug fix: previous code used symmetric T_hop — wrong for p-wave.
            % Fermion antisymmetry demands Δ(i,j) = -Δ(j,i).
            if nargin < 4, bc = 'obc'; end
            sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
            dx = dvec(1);  dy = dvec(2);  dz = dvec(3);

            % d·σ·iσ_y spin matrix
            dsig_isy = [ dz,         dx - 1i*dy ;
                         dx + 1i*dy, -dz         ];

            % Antisymmetric bond matrix A_N
            A_N = zeros(N);
            for i = 1:N-1
                A_N(i,   i+1) = +Delta_t;   % forward bond
                A_N(i+1, i  ) = -Delta_t;   % backward bond (antisymmetric)
            end
            if strcmp(bc,'pbc')
                A_N(N, 1) = +Delta_t;
                A_N(1, N) = -Delta_t;
            end

            D = kron(dsig_isy, A_N);
        end

        % ----------------------------------------------------------------
        % BdG ASSEMBLY
        % ----------------------------------------------------------------
        function H = assemble(H_N, Delta_hat)
            % Assemble 4N×4N BdG Hamiltonian from normal and pairing blocks.
            %
            % H_BdG = [  H_N,     Δ̂   ]
            %          [ -Δ̂†,  -H_N* ]
            %
            % Particle-hole symmetry: Ξ H_BdG Ξ⁻¹ = -H_BdG  ✓
            H = [  H_N,            Delta_hat   ;
                  -Delta_hat',    -conj(H_N)   ];
        end

        % ----------------------------------------------------------------
        % INDEX HELPERS
        % ----------------------------------------------------------------
        function idx = elec_up(N),   idx = 1:N;        end
        function idx = elec_dn(N),   idx = (1:N)+N;    end
        function idx = hole_up(N),   idx = (1:N)+2*N;  end
        function idx = hole_dn(N),   idx = (1:N)+3*N;  end
        function idx = elec_block(N), idx = 1:2*N;     end
        function idx = hole_block(N), idx = (1:2*N)+2*N; end

        % ----------------------------------------------------------------
        % SYMMETRY CHECKS
        % ----------------------------------------------------------------
        function [ok, res] = ph_sym(H)
            % Check particle-hole symmetry: τ_x conj(H) τ_x = -H
            dim = size(H,1);  N2 = dim/2;
            Tx  = [zeros(N2), eye(N2); eye(N2), zeros(N2)];
            res = norm(Tx*conj(H)*Tx + H,'fro') / norm(H,'fro');
            ok  = res < 1e-10;
            if ok
                fprintf('  PH symmetry : OK  (residual=%.2e)\n', res);
            else
                fprintf('  PH symmetry : VIOLATED  (residual=%.2e)\n', res);
            end
        end

        function [ok, res] = tr_sym(H)
            % Check time-reversal symmetry: (iσ_y⊗I) H* (-iσ_y⊗I) = H
            dim = size(H,1);  N4 = dim/4;
            sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
            isy = kron(1i*sig.sy, eye(2*N4));
            res = norm(isy*conj(H)*isy' - H,'fro') / (norm(H,'fro') + eps);
            ok  = res < 1e-10;
            if ok
                fprintf('  TR symmetry : OK  (residual=%.2e)\n', res);
            else
                fprintf('  TR symmetry : broken  (residual=%.2e)\n', res);
            end
        end

        % ----------------------------------------------------------------
        % INFO
        % ----------------------------------------------------------------
        function info()
            fprintf('\n  UQMaT Nambu basis: Ψ = (c↑, c↓, c↑†, c↓†)ᵀ\n');
            fprintf('  ------------------------------------------------\n');
            fprintf('  Block      Operator   Rows (N sites)\n');
            fprintf('  ---------  ---------  -----------------\n');
            fprintf('  Electron↑  c_{i↑}    1      .. N\n');
            fprintf('  Electron↓  c_{i↓}    N+1    .. 2N\n');
            fprintf('  Hole↑      c†_{i↑}   2N+1   .. 3N\n');
            fprintf('  Hole↓      c†_{i↓}   3N+1   .. 4N\n');
            fprintf('\n');
            fprintf('  H_BdG = [  H_N,     Δ̂   ]\n');
            fprintf('          [ -Δ̂†,  -H_N* ]\n\n');
            fprintf('  Singlet: Δ̂ = Δ_s·(iσ_y) ⊗ I_N\n');
            fprintf('  Triplet: Δ̂ = (d⃗·σ⃗)(iσ_y) ⊗ A_N   [A antisymmetric]\n\n');
        end
    end
end
