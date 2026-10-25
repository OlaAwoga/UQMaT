classdef Geometry
    % UQMaT.Models.Lattice2D.Kagome.Geometry
    %   Site indexing and geometry utilities for the 2D kagome lattice.
    %
    %   Site index: n(ix,iy,sub) = 3*((ix-1)*Ny + (iy-1)) + sub
    %               sub: 1=A, 2=B, 3=C
    %
    %   Sublattice positions within unit cell (fractional coords):
    %     A: (0,    0   )
    %     B: (0.5,  0   )
    %     C: (0.25, 0.5 )
    %
    %   BZ path: Γ=(0,0) → M=(π,0) → K=(2π/3, 2π/3√3) → Γ

    methods (Static)

        function n = idx(ix, iy, sub, Ny)
            % Site index (1-based)
            n = 3*((ix-1)*Ny + (iy-1)) + sub;
        end

        function [ix, iy, sub] = inv_idx(n, Ny)
            % Inverse: n → (ix, iy, sub)
            n   = n - 1;
            sub = mod(n, 3) + 1;
            nc  = (n - sub + 1) / 3;
            ix  = floor(nc / Ny) + 1;
            iy  = mod(nc, Ny) + 1;
        end

        function [xyz] = site_positions(Nx, Ny)
            % Return Ns×2 array of (x,y) positions in units of the
            % lattice constant a (triangle side length).
            % Lattice vectors: a1=[2,0], a2=[1,sqrt(3)]
            % Sublattice offsets: A=[0,0], B=[1,0], C=[0.5, sqrt(3)/2]
            Ns  = 3*Nx*Ny;
            xyz = zeros(Ns, 2);
            sub_offset = [0, 0;
                          1, 0;
                          0.5, sqrt(3)/2];
            for ix = 1:Nx
                for iy = 1:Ny
                    r0 = (ix-1)*[2,0] + (iy-1)*[1,sqrt(3)];
                    for sub = 1:3
                        n = UQMaT.Models.Lattice2D.Kagome.Geometry.idx(ix,iy,sub,Ny);
                        xyz(n,:) = r0 + sub_offset(sub,:);
                    end
                end
            end
        end

        function [kpath, klabels, kticks] = kpath_kagome(Nk)
            if nargin < 1, Nk = 300; end
            Ns = ceil(Nk/3);
            G  = [0,       0];
            M  = [pi,      0];
            K  = [2*pi/3,  2*pi/3/sqrt(3)];
            s1 = [linspace(G(1),M(1),Ns)', linspace(G(2),M(2),Ns)'];
            s2 = [linspace(M(1),K(1),Ns)', linspace(M(2),K(2),Ns)'];
            s3 = [linspace(K(1),G(1),Ns)', linspace(K(2),G(2),Ns)'];
            kpath   = [s1; s2(2:end,:); s3(2:end,:)];
            klabels = {'\Gamma','M','K','\Gamma'};
            kticks  = [1, Ns, 2*Ns-1, size(kpath,1)];
        end

        function info(Nx, Ny, bc)
            Ns = 3*Nx*Ny;
            fprintf('\n  Kagome %d×%d=%d sites  bc=%s\n', Nx, Ny, Ns, upper(bc));
            fprintf('  n(ix,iy,sub) = 3*((ix-1)*Ny + (iy-1)) + sub\n');
            fprintf('  sub: 1=A, 2=B, 3=C\n');
            fprintf('  Flat band at E = +2t\n\n');
        end

    end
end
