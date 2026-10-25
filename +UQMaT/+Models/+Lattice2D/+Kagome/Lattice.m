classdef Lattice < UQMaT.Models.Base
    % UQMaT.Models.Lattice2D.Kagome.Lattice
    %   2D kagome tight-binding lattice.
    %
    %   Unit cell: 3 sites (A, B, C) arranged as a corner-sharing triangle.
    %   Total sites: Ns = 3 * Nx * Ny
    %   Site index: n(ix,iy,sub) = 3*((ix-1)*Ny + (iy-1)) + sub
    %               sub: 1=A, 2=B, 3=C
    %
    %   Lattice vectors (in units of the triangle side length a):
    %     a1 = [2, 0]        (x direction)
    %     a2 = [1, sqrt(3)]  (60° from x)
    %
    %   Within unit cell (ix,iy), NN bonds:
    %     A(ix,iy)–B(ix,iy)      intra-cell
    %     A(ix,iy)–C(ix,iy)      intra-cell
    %     B(ix,iy)–C(ix,iy)      intra-cell
    %   Inter-cell NN bonds:
    %     B(ix,iy)–A(ix+1,iy)    (x direction)
    %     C(ix,iy)–A(ix,iy+1)    (y direction)
    %     C(ix,iy)–B(ix-1,iy+1)  (diagonal)
    %
    %   Spectrum: 1 flat band at E=+2t, 2 dispersive bands
    %             Flat band touches upper dispersive band at K,K'
    %
    %   H sizes:
    %     sc=false, spin=false  →  Ns × Ns
    %     sc=false, spin=true   →  2Ns × 2Ns
    %     sc=true               →  4Ns × 4Ns  Ψ=(c↑,c↓,c↑†,c↓†)ᵀ
    %
    %   Usage:
    %     p   = UQMaT.Params('system','kagome','Nx',6,'Ny',6,'t',1,'mu',0)
    %     sys = UQMaT.System(p)
    %     [kd, Ek, kl, kt] = UQMaT.Models.Lattice2D.Kagome.BC.bands(sys)

    %  Citation
    %  --------
    %  Syozi I., Prog. Theor. Phys. 6, 306 (1951). [kagome lattice]
    %  Flat band: Derzhko O. et al., Prog. Theor. Phys. Suppl. 176, 333 (2008).
    %  Topological SC on kagome: Ko W.-H. et al., PRB 79, 214502 (2009).
    %  Chiral SC: Mazin I.I. et al., Nat. Comm. 5, 4261 (2014).
    methods
        function obj = Lattice(p)
            obj.params = p;
            obj.Ndim   = 3 * p.Nx * p.Ny;
            obj.H      = obj.buildH();
        end

        function s = description(obj)
            p  = obj.params;
            Ns = 3*p.Nx*p.Ny;
            if p.sc
                s = sprintf('Kagome SC BdG 4×%d (%d×%d cells) bc=%s', ...
                    Ns, p.Nx, p.Ny, upper(p.bc));
            elseif p.spin
                s = sprintf('Kagome spinful 2×%d (%d×%d cells) bc=%s', ...
                    Ns, p.Nx, p.Ny, upper(p.bc));
            else
                s = sprintf('Kagome normal %d (%d×%d cells) bc=%s', ...
                    Ns, p.Nx, p.Ny, upper(p.bc));
            end
        end

        function H = buildH(obj)
            p   = obj.params;
            Nx  = p.Nx;  Ny = p.Ny;
            Ns  = 3*Nx*Ny;
            bdg = UQMaT.SC.BdG;
            bc  = p.bc;

            % Order parameter contributions
            [D_s_ops, D_t_ops, h_AM, h_CDW] = UQMaT.OrderParameterSymmetry.OPS.build(p);

            % ---- Build kinetic h0 (Ns×Ns) ----
            h0 = obj.hopping_matrix(Nx, Ny, p.t, p.t2, bc);
            h0 = h0 - p.mu*eye(Ns) + p.onsite_val*eye(Ns);

            % Disorder
            h0 = h0 + diag(UQMaT.Disorder.Disorder.kagome_onsite(p, Ns));

            % CDW
            h0 = h0 + h_CDW;

            if p.sc
                H_N = kron(eye(2), h0);
                if any(p.B ~= 0)
                    H_N = H_N + bdg.zeeman(p.B, Ns);
                end
                if p.lambda_R ~= 0
                    A_R = obj.rashba_matrix(Nx, Ny, bc);
                    H_N = H_N + p.lambda_R * kron(bdg.sigma().sy, A_R);
                end
                H_N = H_N + h_AM;
                D   = D_s_ops + D_t_ops;
                H   = bdg.assemble(H_N, D);

            elseif p.spin
                s  = bdg.sigma();
                H  = kron(eye(2), h0);
                if any(p.B ~= 0)
                    H = H + bdg.zeeman(p.B, Ns);
                end
                if p.lambda_R ~= 0
                    A_R = obj.rashba_matrix(Nx, Ny, bc);
                    H   = H + p.lambda_R * kron(s.sy, A_R);
                end
                H = H + h_AM;

            else
                H = h0;
            end
        end
    end

    methods (Static)

        function n = site_index(ix, iy, sub, Ny)
            % n(ix,iy,sub) = 3*((ix-1)*Ny + (iy-1)) + sub
            n = 3*((ix-1)*Ny + (iy-1)) + sub;
        end

        function T = hopping_matrix(Nx, Ny, t, t2, bc)
            % Build Ns×Ns kagome hopping matrix.
            % NN: weight -t  (6 bonds per unit cell)
            % NNN: weight -t2 (if nonzero)
            Ns = 3*Nx*Ny;
            T  = zeros(Ns);

            pbc = strcmp(bc,'pbc');
            idx = @(ix,iy,s) UQMaT.Models.Lattice2D.Kagome.Lattice.idx_wrap(ix,iy,s,Nx,Ny,pbc);

            for ix = 1:Nx
                for iy = 1:Ny
                    nA = idx(ix, iy, 1);
                    nB = idx(ix, iy, 2);
                    nC = idx(ix, iy, 3);

                    % ---- Intra-cell NN bonds (A-B, A-C, B-C) ----
                    if ~isnan(nA) && ~isnan(nB)
                        T(nA,nB) = T(nA,nB) - t;
                        T(nB,nA) = T(nB,nA) - t;
                    end
                    if ~isnan(nA) && ~isnan(nC)
                        T(nA,nC) = T(nA,nC) - t;
                        T(nC,nA) = T(nC,nA) - t;
                    end
                    if ~isnan(nB) && ~isnan(nC)
                        T(nB,nC) = T(nB,nC) - t;
                        T(nC,nB) = T(nC,nB) - t;
                    end

                    % ---- Inter-cell NN bonds ----
                    % B(ix,iy) – A(ix+1,iy)   (x direction)
                    nA2 = idx(ix+1, iy, 1);
                    if ~isnan(nB) && ~isnan(nA2)
                        T(nB,nA2) = T(nB,nA2) - t;
                        T(nA2,nB) = T(nA2,nB) - t;
                    end
                    % C(ix,iy) – A(ix,iy+1)   (y direction)
                    nA3 = idx(ix, iy+1, 1);
                    if ~isnan(nC) && ~isnan(nA3)
                        T(nC,nA3) = T(nC,nA3) - t;
                        T(nA3,nC) = T(nA3,nC) - t;
                    end
                    % C(ix,iy) – B(ix-1,iy+1) (diagonal)
                    nB2 = idx(ix-1, iy+1, 2);
                    if ~isnan(nC) && ~isnan(nB2)
                        T(nC,nB2) = T(nC,nB2) - t;
                        T(nB2,nC) = T(nB2,nC) - t;
                    end

                    % ---- NNN bonds (if t2 ~= 0) ----
                    if t2 ~= 0
                        T = UQMaT.Models.Lattice2D.Kagome.Lattice.add_nnn(...
                            T, ix, iy, Nx, Ny, t2, pbc);
                    end
                end
            end
        end

    end

    methods (Static, Access = private)

        function n = idx_wrap(ix, iy, sub, Nx, Ny, pbc)
            if pbc
                ix = mod(ix-1,Nx)+1;
                iy = mod(iy-1,Ny)+1;
            end
            if ix < 1 || ix > Nx || iy < 1 || iy > Ny
                n = NaN;
            else
                n = UQMaT.Models.Lattice2D.Kagome.Lattice.site_index(ix,iy,sub,Ny);
            end
        end

        function T = add_nnn(T, ix, iy, Nx, Ny, t2, pbc)
            % NNN hoppings: same-sublattice sites two bonds away
            % A-A: (ix±1,iy), (ix,iy±1)
            % B-B: same
            % C-C: same
            idx = @(x,y,s) UQMaT.Models.Lattice2D.Kagome.Lattice.idx_wrap(x,y,s,Nx,Ny,pbc);
            dirs = [+1,0; -1,0; 0,+1; 0,-1];
            for sub = 1:3
                n0 = idx(ix,iy,sub);
                if isnan(n0), continue; end
                for k = 1:4
                    n1 = idx(ix+dirs(k,1), iy+dirs(k,2), sub);
                    if ~isnan(n1)
                        T(n0,n1) = T(n0,n1) - t2;
                    end
                end
            end
        end

        function A = rashba_matrix(Nx, Ny, bc)
            % Antisymmetric NN bond matrix for Rashba SOC on kagome
            Ns  = 3*Nx*Ny;
            A   = zeros(Ns);
            pbc = strcmp(bc,'pbc');
            idx = @(ix,iy,s) UQMaT.Models.Lattice2D.Kagome.Lattice.idx_wrap(ix,iy,s,Nx,Ny,pbc);

            for ix = 1:Nx
                for iy = 1:Ny
                    pairs = {1,2; 1,3; 2,3};  % A-B, A-C, B-C (intra)
                    for k = 1:3
                        n1 = idx(ix,iy,pairs{k,1});
                        n2 = idx(ix,iy,pairs{k,2});
                        if ~isnan(n1) && ~isnan(n2)
                            A(n1,n2) = +1;
                            A(n2,n1) = -1;
                        end
                    end
                    % Inter-cell antisymmetric bonds
                    nB  = idx(ix,iy,2); nA2 = idx(ix+1,iy,1);
                    if ~isnan(nB) && ~isnan(nA2)
                        A(nB,nA2) = +1; A(nA2,nB) = -1;
                    end
                    nC  = idx(ix,iy,3); nA3 = idx(ix,iy+1,1);
                    if ~isnan(nC) && ~isnan(nA3)
                        A(nC,nA3) = +1; A(nA3,nC) = -1;
                    end
                    nC2 = idx(ix,iy,3); nB2 = idx(ix-1,iy+1,2);
                    if ~isnan(nC2) && ~isnan(nB2)
                        A(nC2,nB2) = +1; A(nB2,nC2) = -1;
                    end
                end
            end
        end

    end
end
