classdef BC
    % UQMaT.Models.Lattice.BoundaryConditions.BC
    %   Boundary condition utilities and band structure via twisted BC.
    %
    % ===================================================================
    %  THREE BOUNDARY CONDITIONS
    % ===================================================================
    %
    %   OBC  (open, bc='obc'):
    %     T(1,N) = T(N,1) = 0           free chain ends
    %     Physical: edge states, finite-size effects
    %     k is NOT a good quantum number
    %
    %   PBC  (periodic, bc='pbc'):
    %     T(1,N) = T(N,1) = -t          hard wrap-around bond
    %     Physical: bulk periodic, discrete k_n = 2πn/N
    %
    %   TBC  (twisted / Bloch, bc='tbc', k is a parameter):
    %     T(k)(1,N) = -t · exp(-ik)     phase-twisted wrap bond
    %     T(k)(N,1) = -t · exp(+ik)
    %     Physical: Bloch Hamiltonian — scanning k ∈ [-π,π] gives bands
    %
    % ===================================================================
    %  BAND STRUCTURE VIA TWISTED BOUNDARY CONDITIONS
    % ===================================================================
    %
    %   The Bloch Hamiltonian H(k) is obtained from the OBC Hamiltonian
    %   by replacing the missing wrap-around bond with a phase factor:
    %
    %     H(k) = H_OBC + (-t)·exp(-ik)·|1><N|  +  (-t)·exp(+ik)·|N><1|
    %
    %   This is equivalent to the Fourier transform of the real-space H
    %   in the thermodynamic limit, and does NOT require writing H(k)
    %   separately — it is derived directly from the same H that gives
    %   OBC and PBC.
    %
    %   For a uniform chain:       E(k) = -2t cos(k) - μ     (exact)
    %   For a 2-site SSH unit cell: E±(k) = ±|v + w·exp(-ik)| (exact)
    %   For Kitaev BdG:            as 2×2 Nambu H(k)
    %   For spinful + Rashba:      as 2×2 spin H(k)
    %   For spinful BdG:           as 4×4 Nambu⊗spin H(k)
    %
    %   In 2D, the phase is applied separately to x- and y-bonds:
    %     T_x(kx)(1,Ny): bond, Nx → 1 at fixed iy:  phase exp(-ikx)
    %     T_y(ky)(1,Nx): bond, Ny → 1 at fixed ix:  phase exp(-iky)
    %
    % ===================================================================
    %  USAGE
    % ===================================================================
    %
    %   --- 1D band structure ---
    %   p   = UQMaT.Params('system','chain','N',1,'bc','obc','t',1,'mu',0);
    %   sys = UQMaT.System(p);
    %   [kn, Ek] = UQMaT.Models.Lattice.BoundaryConditions.BC.bands1D(sys);
    %   [kn, Ek] = UQMaT.Models.Lattice.BoundaryConditions.BC.bands1D(sys, 200);
    %
    %   Works for any system size; even N=1 gives the correct dispersion
    %   because H(k) is built from the parameters, not from a large matrix.
    %
    %   --- 2D band structure ---
    %   p   = UQMaT.Params('system','square','Nx',1,'Ny',1,'t',1,'mu',0);
    %   sys = UQMaT.System(p);
    %   [kd, Ek, kl, kt] = UQMaT.Models.Lattice.BoundaryConditions.BC.bands2D(sys);

    methods (Static)

        % ================================================================
        %  1D HOPPING MATRIX  (real-space, OBC or PBC)
        % ================================================================
        function T = hop1D(N, t, t2, bc)
            % Build N×N 1D hopping matrix.
            %   bc = 'obc' : no wrap-around bond  (T(1,N) = 0)
            %   bc = 'pbc' : real wrap-around bond (T(1,N) = -t)
            %   bc = 'tbc' : call hop1D_bloch(N,t,t2,k) instead
            if nargin < 3 || isempty(t2), t2 = 0;     end
            if nargin < 4,                bc = 'obc';  end
            bc = lower(strtrim(bc));

            e1 = -t * ones(N-1,1);
            T  = diag(e1,1) + diag(e1,-1);
            if t2 ~= 0
                e2 = -t2 * ones(N-2,1);
                T  = T + diag(e2,2) + diag(e2,-2);
            end
            if strcmp(bc,'pbc')
                T(1,N) = -t;   T(N,1) = -t;
                if t2 ~= 0 && N > 3
                    T(1,N-1)=-t2; T(2,N)=-t2;
                    T(N-1,1)=-t2; T(N,2)=-t2;
                end
            end
            % 'tbc': return plain OBC matrix; caller adds phase in hop1D_bloch
        end

        % ================================================================
        %  1D BLOCH HOPPING MATRIX  (twisted BC at wavevector k)
        % ================================================================
        function T = hop1D_bloch(N, t, t2, k)
            % Build N×N 1D hopping matrix with twisted (Bloch) BC.
            %
            %   T(k)(i,i+1) = -t          (bulk bonds, unchanged)
            %   T(k)(1, N)  = -t·exp(-ik) (wrap bond, forward, gains phase)
            %   T(k)(N, 1)  = -t·exp(+ik) (wrap bond, backward, loses phase)
            %
            % The phase exp(ik) comes from translating by one lattice constant.
            % Scanning k ∈ [-π,π] and diagonalising T(k) gives the bulk bands.
            if nargin < 3 || isempty(t2), t2 = 0; end

            % Start from OBC matrix
            bc = UQMaT.Models.Lattice.BoundaryConditions.BC;
            T  = bc.hop1D(N, t, t2, 'obc');

            % Add twisted wrap-around bond
            T(1, N) = -t * exp(-1i*k);
            T(N, 1) = -t * exp(+1i*k);

            % Next-nearest twisted bonds (second wrap-around)
            if t2 ~= 0 && N > 2
                T(1, N-1) = -t2 * exp(-1i*k);
                T(N-1, 1) = -t2 * exp(+1i*k);
                T(2, N)   = -t2 * exp(-1i*k);
                T(N, 2)   = -t2 * exp(+1i*k);
            end
        end

        % ================================================================
        %  2D HOPPING MATRIX  (real-space, OBC or PBC)
        % ================================================================
        function T = hop2D(Nx, Ny, t, t2, bcx, bcy)
            % Build Ns×Ns 2D square lattice hopping matrix (sparse).
            % T = kron(Tx, I_y) + kron(I_x, Ty)
            if nargin < 4 || isempty(t2), t2  = 0;     end
            if nargin < 5,                bcx = 'obc';  end
            if nargin < 6,                bcy = 'obc';  end
            bc = UQMaT.Models.Lattice.BoundaryConditions.BC;
            Tx = bc.hop1D(Nx, t, t2, bcx);
            Ty = bc.hop1D(Ny, t, t2, bcy);
            T  = kron(Tx, speye(Ny)) + kron(speye(Nx), Ty);
        end

        % ================================================================
        %  2D BLOCH HOPPING MATRIX  (twisted BC at wavevector kx, ky)
        % ================================================================
        function T = hop2D_bloch(Nx, Ny, t, t2, kx, ky)
            % Build Ns×Ns 2D hopping matrix with twisted BC.
            %
            % x-bonds: T_x(kx) — wrap bond at ix=Nx→ix=1 gains phase exp(-ikx)
            % y-bonds: T_y(ky) — wrap bond at iy=Ny→iy=1 gains phase exp(-iky)
            % Then: T = kron(Tx_bloch, I_y) + kron(I_x, Ty_bloch)
            if nargin < 4 || isempty(t2), t2 = 0; end
            bc = UQMaT.Models.Lattice.BoundaryConditions.BC;
            Tx = bc.hop1D_bloch(Nx, t, t2, kx);
            Ty = bc.hop1D_bloch(Ny, t, t2, ky);
            T  = kron(Tx, speye(Ny)) + kron(speye(Nx), Ty);
        end

        % ================================================================
        %  BAND STRUCTURE  1D
        % ================================================================
        function [kn, Ek] = bands1D(sys, Nk)
            % Extract 1D band structure via twisted boundary conditions.
            %
            %   For each k in [-π, π], build H(k) by replacing the
            %   wrap-around bond of H_OBC with -t·exp(±ik), then diagonalise.
            %   No separate H(k) expression is written — the parameters
            %   are read directly from sys.params.
            %
            %   [kn, Ek] = bands1D(sys)        % 200 k-points (default)
            %   [kn, Ek] = bands1D(sys, Nk)    % Nk k-points
            %
            %   kn : Nk×1   k-values in [-π, π]
            %   Ek : Nk×Nb  band energies  (Nb = 1 normal, 2 spin, 4 BdG)
            %
            %   Works for: chain, ssh, kitaev (all sizes including N=1 unit cell)
            if nargin < 2, Nk = 200; end
            p  = sys.params;
            bc = UQMaT.Models.Lattice.BoundaryConditions.BC;
            kn = linspace(-pi, pi, Nk)';

            switch p.model_class
                case 'regular'
                    Ek = bc.bands1D_regular(p, kn);
                case 'ssh'
                    Ek = bc.bands1D_ssh(p, kn);
                case 'kitaev'
                    Ek = bc.bands1D_kitaev(p, kn);
                otherwise
                    error('UQMaT:BC:bands1D:unsupported', ...
                        'bands1D not supported for model_class="%s".', ...
                        p.model_class);
            end
        end

        % ================================================================
        %  BAND STRUCTURE  2D  (along Γ→X→M→Γ)
        % ================================================================
        function [kdist, Ek, klabels, kticks] = bands2D(sys, Nk)
            % Extract 2D band structure via twisted boundary conditions.
            %
            % The k-path and meaning depend on which directions are periodic:
            %
            %   bcx='pbc', bcy='pbc'  →  full 2D path  Γ→X→M→Γ
            %                             both kx and ky are good QNs
            %                             kpts is Nk × 2
            %
            %   bcx='pbc', bcy='obc'  →  1D path  kx ∈ [-π,π], ky=0 (fixed)
            %                             only kx is a good QN
            %                             cylinder geometry (periodic in x)
            %
            %   bcx='obc', bcy='pbc'  →  1D path  ky ∈ [-π,π], kx=0 (fixed)
            %                             only ky is a good QN
            %                             cylinder geometry (periodic in y)
            %
            %   bcx='obc', bcy='obc'  →  error: no good k in either direction.
            %                             Use sys.E for the real-space spectrum.
            %
            % For mixed BC (one periodic, one open), the H(k) is built by
            % twisting only the periodic wrap bond; the open direction retains
            % its full Ny (or Nx) real-space structure. This gives bands
            % dispersing along the periodic direction with the confined-
            % direction sub-band structure already encoded.
            %
            %   [kdist, Ek, klabels, kticks] = bands2D(sys)       % 200 k-pts
            %   [kdist, Ek, klabels, kticks] = bands2D(sys, Nk)

            if nargin < 2, Nk = 200; end
            p  = sys.params;
            bc = UQMaT.Models.Lattice.BoundaryConditions.BC;

            pbc_x = strcmp(p.bcx, 'pbc');
            pbc_y = strcmp(p.bcy, 'pbc');

            if pbc_x && pbc_y
                % ---- Full 2D: Γ→X→M→Γ -----------------------------------
                Ns    = ceil(Nk/3);
                G     = [0,0];  X = [pi,0];  M = [pi,pi];
                seg1  = [linspace(G(1),X(1),Ns)', linspace(G(2),X(2),Ns)'];
                seg2  = [linspace(X(1),M(1),Ns)', linspace(X(2),M(2),Ns)'];
                seg3  = [linspace(M(1),G(1),Ns)', linspace(M(2),G(2),Ns)'];
                kpts  = [seg1; seg2(2:end,:); seg3(2:end,:)];
                klabels = {'\Gamma', 'X', 'M', '\Gamma'};
                kticks  = [1, Ns, 2*Ns-1, size(kpts,1)];
                Ek      = bc.bands2D_square(p, kpts, 'xy');

            elseif pbc_x && ~pbc_y
                % ---- Cylinder x: kx ∈ [-π,π], ky fixed = 0 --------------
                % H(kx) is Ny×Ny (or Nb*Ny) in the open y-direction
                kx_vec  = linspace(-pi, pi, Nk)';
                kpts    = [kx_vec, zeros(Nk,1)];
                klabels = {'-\pi', '0', '\pi'};
                kticks  = [1, round(Nk/2), Nk];
                Ek      = bc.bands2D_square(p, kpts, 'x');

            elseif ~pbc_x && pbc_y
                % ---- Cylinder y: ky ∈ [-π,π], kx fixed = 0 --------------
                ky_vec  = linspace(-pi, pi, Nk)';
                kpts    = [zeros(Nk,1), ky_vec];
                klabels = {'-\pi', '0', '\pi'};
                kticks  = [1, round(Nk/2), Nk];
                Ek      = bc.bands2D_square(p, kpts, 'y');

            else
                error('UQMaT:BC:bands2D:noGoodK', ...
                    ['bcx=''obc'' and bcy=''obc'': no good k in either direction.\n' ...
                     'Use sys.E to see the real-space energy spectrum.']);
            end

            % Cumulative k-path distance
            kdist = zeros(size(kpts,1), 1);
            for ik = 2:size(kpts,1)
                kdist(ik) = kdist(ik-1) + norm(kpts(ik,:)-kpts(ik-1,:));
            end
        end

        % ================================================================
        %  INFO
        % ================================================================
        function info()
            fprintf('\n  UQMaT Boundary Conditions\n');
            fprintf('  ==========================\n\n');
            fprintf('  OBC  (bc=''obc'')  : T(1,N) = T(N,1) = 0\n');
            fprintf('                      free chain ends\n');
            fprintf('                      use for: edge states, finite chains\n\n');
            fprintf('  PBC  (bc=''pbc'')  : T(1,N) = T(N,1) = -t\n');
            fprintf('                      hard periodic bond\n');
            fprintf('                      use for: bulk, avoid surface effects\n\n');
            fprintf('  TBC  (Bloch)      : T(k)(1,N) = -t·exp(-ik)\n');
            fprintf('                      T(k)(N,1) = -t·exp(+ik)\n');
            fprintf('                      scan k ∈ [-π,π] → band structure\n');
            fprintf('                      same H as OBC + phase on wrap bond\n\n');
            fprintf('  bands1D(sys) / bands2D(sys) compute E(k) this way.\n\n');
        end

    end   % public static methods

    % ====================================================================
    %  Private band structure builders (one per model class)
    % ====================================================================
    methods (Static, Access = private)

        % ---- Regular chain (normal / spin / BdG) -----------------------
        function Ek = bands1D_regular(p, kn)
            bc  = UQMaT.Models.Lattice.BoundaryConditions.BC;
            nb  = UQMaT.Models.Lattice.Nambu.Basis;
            sig = nb.sigma();
            Nk  = numel(kn);

            % Determine number of bands
            if p.bdg,       Nb = 4;
            elseif p.spin,  Nb = 2;
            else,           Nb = 1;
            end
            Ek = zeros(Nk, Nb);

            for ik = 1:Nk
                k = kn(ik);

                % Bloch hopping (1-site unit cell → scalar ε(k))
                % T_bloch(k) for N=1: only the wrap bond exists
                %   ε(k) = -t(e^ik + e^-ik) - t2(e^2ik + e^-2ik) - mu
                %         = -2t cos(k) - 2t2 cos(2k) - mu
                ek = -2*p.t*cos(k) - 2*p.t2*cos(2*k) - p.mu;

                if ~p.spin && ~p.bdg
                    Ek(ik,1) = ek;

                elseif p.spin && ~p.bdg
                    % 2×2 spin Hamiltonian
                    HkN = ek*eye(2) ...
                        + p.B(3)*sig.sz + p.B(1)*sig.sx + p.B(2)*sig.sy;
                    if p.lambda_R ~= 0
                        % Rashba: H_R(k) = λ_R(-sin(k)·σ_y)  (1D x-direction)
                        HkN = HkN - p.lambda_R * sin(k) * sig.sy;
                    end
                    Ek(ik,:) = sort(real(eig(HkN)))';

                else
                    % 4×4 BdG Nambu Hamiltonian
                    HkN = ek*eye(2) ...
                        + p.B(3)*sig.sz + p.B(1)*sig.sx + p.B(2)*sig.sy;
                    if p.lambda_R ~= 0
                        HkN = HkN - p.lambda_R * sin(k) * sig.sy;
                    end
                    % Pairing block Δ̂(k)
                    Dk = bc.pairing_bloch(p, k);
                    % BdG: [HkN, Dk; -Dk†, -HkN*]
                    Hbdg = nb.assemble(HkN, Dk);
                    Ek(ik,:) = sort(real(eig(Hbdg)))';
                end
            end
        end

        % ---- SSH chain ------------------------------------------------
        function Ek = bands1D_ssh(p, kn)
            % SSH: 2-site unit cell (A, B)
            % H(k) = [mu,  v + w·exp(-ik); v + w·exp(+ik),  mu]
            Nk = numel(kn);
            Ek = zeros(Nk, 2);
            for ik = 1:Nk
                k = kn(ik);
                off = p.v + p.w * exp(-1i*k);
                Hk  = [p.mu, off; conj(off), p.mu];
                Ek(ik,:) = sort(real(eig(Hk)))';
            end
        end

        % ---- Kitaev chain (spinless p-wave BdG) -----------------------
        function Ek = bands1D_kitaev(p, kn)
            % Kitaev H(k) — 2×2 Nambu, spinless
            % H(k) = [-2t cos(k) - mu,   2i·Δ_p·sin(k) ]
            %        [-2i·Δ_p·sin(k),    2t cos(k) + mu ]
            Nk  = numel(kn);
            Dp  = mean(p.sc_Delta_s_p_vec);
            if Dp == 0
                % Normal state (no BdG)
                Ek = zeros(Nk,1);
                for ik=1:Nk
                    Ek(ik) = -2*p.t*cos(kn(ik)) - p.mu;
                end
            else
                Ek = zeros(Nk, 2);
                for ik = 1:Nk
                    k  = kn(ik);
                    ek = -2*p.t*cos(k) - p.mu;
                    Hk = [ek,             2i*Dp*sin(k);
                          -2i*Dp*sin(k), -ek           ];
                    Ek(ik,:) = sort(real(eig(Hk)))';
                end
            end
        end

        % ---- 2D square lattice (normal / spin / BdG) ------------------
        function Ek = bands2D_square(p, kpts, dir)
            % dir = 'xy' : both directions periodic → 1x1 Bloch H(kx,ky)
            %       'x'  : only x periodic → Ny×Ny H(kx) with full y structure
            %       'y'  : only y periodic → Nx×Nx H(ky) with full x structure
            bc  = UQMaT.Models.Lattice.BoundaryConditions.BC;
            nb  = UQMaT.Models.Lattice.Nambu.Basis;
            sig = nb.sigma();
            Nkpts = size(kpts, 1);

            if nargin < 3, dir = 'xy'; end

            switch dir

                % ============================================================
                case 'xy'
                % Both periodic: scalar ε(kx,ky), small Bloch H
                % ============================================================
                    if p.bdg,      Nb = 4;
                    elseif p.spin, Nb = 2;
                    else,          Nb = 1;
                    end
                    Ek = zeros(Nkpts, Nb);

                    for ik = 1:Nkpts
                        kx = kpts(ik,1);  ky = kpts(ik,2);
                        ek = -2*p.t*(cos(kx)+cos(ky)) ...
                             - 4*p.t2*cos(kx)*cos(ky) - p.mu;

                        if ~p.spin && ~p.bdg
                            Ek(ik,1) = ek;
                        elseif p.spin && ~p.bdg
                            HkN = ek*eye(2) ...
                                + p.B(3)*sig.sz + p.B(1)*sig.sx + p.B(2)*sig.sy;
                            if p.lambda_R ~= 0
                                HkN = HkN + p.lambda_R * ...
                                    (sin(ky)*sig.sx - sin(kx)*sig.sy);
                            end
                            Ek(ik,:) = sort(real(eig(HkN)))';
                        else
                            HkN = ek*eye(2) ...
                                + p.B(3)*sig.sz + p.B(1)*sig.sx + p.B(2)*sig.sy;
                            if p.lambda_R ~= 0
                                HkN = HkN + p.lambda_R * ...
                                    (sin(ky)*sig.sx - sin(kx)*sig.sy);
                            end
                            Dk   = bc.pairing_bloch(p, kx, ky);
                            Hbdg = nb.assemble(HkN, Dk);
                            Ek(ik,:) = sort(real(eig(Hbdg)))';
                        end
                    end

                % ============================================================
                case 'x'
                % x periodic, y open: H(kx) is Ny×Ny (open y chain of Ny sites)
                % Each kx gives a different Ny×Ny normal Hamiltonian:
                %   H(kx)_OBC_y = -2t cos(kx)·I_Ny + T_y_OBC - mu·I_Ny
                %   (the kx contribution shifts the on-site energy uniformly)
                % Rashba couples the spin blocks across y-bonds only
                % ============================================================
                    Ny   = p.Ny;
                    Nblk = bc.nblk(p);
                    Ek   = zeros(Nkpts, Nblk * Ny);

                    % y hopping (OBC, real-space Ny×Ny)
                    Ty_obc = bc.hop1D(Ny, p.t, 0, 'obc');
                    if p.lambda_R ~= 0
                        Ry_obc = bc.hop1D(Ny, p.lambda_R, 0, 'obc');
                    end

                    for ik = 1:Nkpts
                        kx = kpts(ik,1);

                        % On-site energy from kx dispersion + chemical potential
                        eps_kx = -2*p.t*cos(kx) - p.mu;

                        % Normal block H_N (Nblk_sp * Ny × Nblk_sp * Ny)
                        if ~p.spin && ~p.bdg
                            H_kx = eps_kx*eye(Ny) + Ty_obc;
                        else
                            % Spin-diagonal hopping in y, kx-shifted on-site
                            H_N = kron(eye(2), eps_kx*eye(Ny) + Ty_obc);

                            % Zeeman (no k-dependence)
                            if p.B(1)~=0||p.B(2)~=0||p.B(3)~=0
                                H_N = H_N + kron( ...
                                    p.B(3)*sig.sz+p.B(1)*sig.sx+p.B(2)*sig.sy, eye(Ny));
                            end

                            % Rashba: x-bond contributes -λ sin(kx) σ_y (scalar)
                            %         y-bond contributes +λ σ_x T_y (real-space)
                            if p.lambda_R ~= 0
                                H_N = H_N ...
                                    - p.lambda_R*sin(kx)*kron(sig.sy, eye(Ny)) ...
                                    + p.lambda_R*kron(sig.sx, Ry_obc);
                            end

                            if ~p.bdg
                                H_kx = H_N;
                            else
                                % Pairing: singlet is on-site → k-indep Ny×Ny block
                                %          triplet has sin(kx) from x-bond
                                Ds = mean(p.sc_Delta_s_s_vec);
                                Dk_s = Ds * kron(1i*sig.sy, eye(Ny));
                                if ismember(p.pairing_mode_mode,{'triplet','mixed'})
                                    dx=[0,0,0](1); dy=[0,0,0](2); dz=[0,0,0](3);
                                    dsisy = [dz,dx-1i*dy; dx+1i*dy,-dz];
                                    % x-bond triplet: sin(kx) factor
                                    Dk_t = p.sc_Delta_s_t_val * sin(kx) * ...
                                           kron(dsisy, eye(Ny));
                                else
                                    Dk_t = zeros(2*Ny);
                                end
                                switch p.pairing_mode_mode
                                    case 'singlet', Dhat = Dk_s;
                                    case 'triplet', Dhat = Dk_t;
                                    case 'mixed',   Dhat = Dk_s + Dk_t;
                                end
                                H_kx = nb.assemble(H_N, Dhat);
                            end
                        end

                        ev = sort(real(eig(H_kx)))';
                        Ek(ik,:) = ev;
                    end

                % ============================================================
                case 'y'
                % y periodic, x open: H(ky) is Nx×Nx (open x chain of Nx sites)
                % Mirror image of the 'x' case with kx ↔ ky, Nx ↔ Ny
                % ============================================================
                    Nx   = p.Nx;
                    Nblk = bc.nblk(p);
                    Ek   = zeros(Nkpts, Nblk * Nx);

                    Tx_obc = bc.hop1D(Nx, p.t, 0, 'obc');
                    if p.lambda_R ~= 0
                        Rx_obc = bc.hop1D(Nx, p.lambda_R, 0, 'obc');
                    end

                    for ik = 1:Nkpts
                        ky = kpts(ik,2);

                        eps_ky = -2*p.t*cos(ky) - p.mu;

                        if ~p.spin && ~p.bdg
                            H_ky = eps_ky*eye(Nx) + Tx_obc;
                        else
                            H_N = kron(eye(2), eps_ky*eye(Nx) + Tx_obc);

                            if p.B(1)~=0||p.B(2)~=0||p.B(3)~=0
                                H_N = H_N + kron( ...
                                    p.B(3)*sig.sz+p.B(1)*sig.sx+p.B(2)*sig.sy, eye(Nx));
                            end

                            if p.lambda_R ~= 0
                                % y-bond: sin(ky) σ_x contribution (scalar)
                                % x-bond: -σ_y T_x (real-space)
                                H_N = H_N ...
                                    + p.lambda_R*sin(ky)*kron(sig.sx, eye(Nx)) ...
                                    - p.lambda_R*kron(sig.sy, Rx_obc);
                            end

                            if ~p.bdg
                                H_ky = H_N;
                            else
                                Ds = mean(p.sc_Delta_s_s_vec);
                                Dk_s = Ds * kron(1i*sig.sy, eye(Nx));
                                if ismember(p.pairing_mode_mode,{'triplet','mixed'})
                                    dx=[0,0,0](1); dy=[0,0,0](2); dz=[0,0,0](3);
                                    dsisy = [dz,dx-1i*dy; dx+1i*dy,-dz];
                                    Dk_t = p.sc_Delta_s_t_val * sin(ky) * ...
                                           kron(dsisy, eye(Nx));
                                else
                                    Dk_t = zeros(2*Nx);
                                end
                                switch p.pairing_mode_mode
                                    case 'singlet', Dhat = Dk_s;
                                    case 'triplet', Dhat = Dk_t;
                                    case 'mixed',   Dhat = Dk_s + Dk_t;
                                end
                                H_ky = nb.assemble(H_N, Dhat);
                            end
                        end

                        ev = sort(real(eig(H_ky)))';
                        Ek(ik,:) = ev;
                    end

            end  % switch dir
        end

        % ---- Helper: number of spin/Nambu blocks -----------------------
        function n = nblk(p)
            if p.bdg,       n = 4;
            elseif p.spin,  n = 2;
            else,           n = 1;
            end
        end

        % ---- Pairing block Δ̂(k)  (2×2 in spin space) ----------------
        function Dk = pairing_bloch(p, kx, ky)
            % Build the 2×2 pairing matrix at wavevector (kx[,ky]).
            %
            %   Singlet s-wave:  Δ̂(k) = Δ_s · iσ_y   (k-independent)
            %   Triplet p-wave:  Δ̂(k) = (d⃗·σ⃗)(iσ_y) · f(k)
            %     1D: f(k)   = sin(k)      (p-wave: odd under k→-k)
            %     2D: fx(kx) = sin(kx), fy(ky) = sin(ky)
            %
            % This is the Bloch-space version of the real-space antisymmetric
            % bond matrix A_N: Fourier transform of A_N(i,i+1)=+Δ,-Δ gives sin(k).
            if nargin < 3, ky = 0; end
            nb  = UQMaT.Models.Lattice.Nambu.Basis;
            sig = nb.sigma();

            switch p.pairing_mode_mode
                case 'singlet'
                    Ds  = mean(p.sc_Delta_s_s_vec);
                    Dk  = Ds * (1i * sig.sy);

                case 'triplet'
                    dx = [0,0,0](1); dy = [0,0,0](2); dz = [0,0,0](3);
                    % d·σ·iσ_y
                    dsig_isy = [dz,         dx - 1i*dy;
                                dx + 1i*dy, -dz        ];
                    % p-wave form factor: sin(k) in 1D; sin(kx)+i·sin(ky) in 2D
                    if ky == 0
                        fk = sin(kx);
                    else
                        fk = sin(kx) + 1i*sin(ky);
                    end
                    Dk = p.sc_Delta_s_t_val * dsig_isy * fk;

                case 'mixed'
                    Ds  = mean(p.sc_Delta_s_s_vec);
                    D_s = Ds * (1i * sig.sy);
                    dx = [0,0,0](1); dy = [0,0,0](2); dz = [0,0,0](3);
                    dsig_isy = [dz,         dx - 1i*dy;
                                dx + 1i*dy, -dz        ];
                    if ky == 0
                        fk = sin(kx);
                    else
                        fk = sin(kx) + 1i*sin(ky);
                    end
                    D_t = p.sc_Delta_s_t_val * dsig_isy * fk;
                    Dk  = D_s + D_t;

                otherwise
                    Dk = zeros(2);
            end
        end

    end   % private static methods
end
