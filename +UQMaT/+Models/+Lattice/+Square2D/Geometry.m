classdef Geometry
    % UQMaT.Models.Lattice.Square2D.Geometry
    %   Index conventions and geometry utilities for the 2D square lattice.
    %
    % ===================================================================
    %  SITE INDEXING
    % ===================================================================
    %
    %   Lattice:  Nx columns (x-direction) × Ny rows (y-direction)
    %   Total sites: Ns = Nx * Ny
    %
    %   Site index (1-based, column-major):
    %     n(ix, iy) = (ix-1)*Ny + iy
    %
    %     ix = 1..Nx   (x-coordinate, column)
    %     iy = 1..Ny   (y-coordinate, row)
    %
    %   Lattice vectors:
    %     a1 = (1, 0)   (x-direction)
    %     a2 = (0, 1)   (y-direction)
    %
    %   Nearest neighbours of site (ix, iy):
    %     (ix±1, iy)   x-bonds   (Nx-1)*Ny bonds + Ny PBC bonds
    %     (ix, iy±1)   y-bonds   Nx*(Ny-1) bonds + Nx PBC bonds
    %
    %   Next-nearest (diagonal) neighbours:
    %     (ix±1, iy±1)   4 bonds per site
    %
    % ===================================================================
    %  NAMBU BLOCK INDICES  (for 4-component BdG, 4Ns × 4Ns)
    % ===================================================================
    %
    %   Ψ = (c↑, c↓, c↑†, c↓†)ᵀ  per site
    %
    %   rows  1    : Ns    →  c_{n↑}    (electron spin-up)
    %   rows  Ns+1 : 2Ns   →  c_{n↓}    (electron spin-down)
    %   rows  2Ns+1: 3Ns   →  c†_{n↑}  (hole spin-up)
    %   rows  3Ns+1: 4Ns   →  c†_{n↓}  (hole spin-down)
    %
    % ===================================================================
    %  HIGH-SYMMETRY K-POINTS  (square BZ)
    % ===================================================================
    %
    %   Γ = (0,  0)
    %   X = (π,  0)
    %   M = (π,  π)
    %   Y = (0,  π)
    %
    %   Standard path: Γ → X → M → Γ

    methods (Static)

        % ---- Site index (ix,iy) → linear index n --------------------
        function n = idx(ix, iy, Ny)
            n = (ix-1)*Ny + iy;
        end

        % ---- Linear index n → (ix, iy) ------------------------------
        function [ix, iy] = inv_idx(n, Ny)
            ix = ceil(n / Ny);
            iy = n - (ix-1)*Ny;
        end

        % ---- All (ix,iy) pairs as Ns×2 array ------------------------
        function coords = all_coords(Nx, Ny)
            Ns     = Nx * Ny;
            coords = zeros(Ns, 2);
            for n = 1:Ns
                [coords(n,1), coords(n,2)] = ...
                    UQMaT.Models.Lattice.Square2D.Geometry.inv_idx(n, Ny);
            end
        end

        % ---- Nearest-neighbour pairs (bonds) ------------------------
        function bonds = nn_bonds(Nx, Ny, bcx, bcy)
            % Returns M×2 array of [n_i, n_j] nearest-neighbour pairs.
            idx = @(ix,iy) UQMaT.Models.Lattice.Square2D.Geometry.idx(ix,iy,Ny);
            bonds = [];
            for ix = 1:Nx
                for iy = 1:Ny
                    n = idx(ix,iy);
                    % x-bond: (ix,iy) → (ix+1,iy)
                    if ix < Nx
                        bonds(end+1,:) = [n, idx(ix+1,iy)]; %#ok<AGROW>
                    elseif strcmp(bcx,'pbc')
                        bonds(end+1,:) = [n, idx(1,iy)];    %#ok<AGROW>
                    end
                    % y-bond: (ix,iy) → (ix,iy+1)
                    if iy < Ny
                        bonds(end+1,:) = [n, idx(ix,iy+1)]; %#ok<AGROW>
                    elseif strcmp(bcy,'pbc')
                        bonds(end+1,:) = [n, idx(ix,1)];    %#ok<AGROW>
                    end
                end
            end
        end

        % ---- Next-nearest-neighbour pairs (diagonal bonds) ----------
        function bonds = nnn_bonds(Nx, Ny, bcx, bcy)
            idx   = @(ix,iy) UQMaT.Models.Lattice.Square2D.Geometry.idx(ix,iy,Ny);
            wrapx = @(ix) mod(ix-1,Nx)+1;
            wrapy = @(iy) mod(iy-1,Ny)+1;
            bonds = [];
            for ix = 1:Nx
                for iy = 1:Ny
                    n = idx(ix,iy);
                    % (+1,+1)
                    jx=ix+1; jy=iy+1;
                    if (jx<=Nx||strcmp(bcx,'pbc')) && (jy<=Ny||strcmp(bcy,'pbc'))
                        bonds(end+1,:)=[n, idx(wrapx(jx),wrapy(jy))]; %#ok<AGROW>
                    end
                    % (+1,-1)
                    jx=ix+1; jy=iy-1;
                    if (jx<=Nx||strcmp(bcx,'pbc')) && (jy>=1||strcmp(bcy,'pbc'))
                        bonds(end+1,:)=[n, idx(wrapx(jx),wrapy(jy))]; %#ok<AGROW>
                    end
                end
            end
        end

        % ---- High-symmetry k-path for band structure ----------------
        function [kpath, klabels, kticks] = kpath_square(Nk)
            % Returns Nk-point path Γ→X→M→Γ through square BZ.
            if nargin < 1, Nk = 100; end
            Ns = ceil(Nk/3);

            G = [0,  0];
            X = [pi, 0];
            M = [pi, pi];

            seg1 = [linspace(G(1),X(1),Ns)', linspace(G(2),X(2),Ns)'];
            seg2 = [linspace(X(1),M(1),Ns)', linspace(X(2),M(2),Ns)'];
            seg3 = [linspace(M(1),G(1),Ns)', linspace(M(2),G(2),Ns)'];

            kpath   = [seg1; seg2(2:end,:); seg3(2:end,:)];
            klabels = {'\Gamma','X','M','\Gamma'};
            kticks  = [1, Ns, 2*Ns-1, size(kpath,1)];
        end

        % ---- Visualise lattice (requires figure) --------------------
        function plot_lattice(Nx, Ny, highlight_sites)
            % Draw the lattice sites and bonds.
            if nargin < 3, highlight_sites = []; end
            coords = UQMaT.Models.Lattice.Square2D.Geometry.all_coords(Nx,Ny);
            hold on;
            % Sites
            plot(coords(:,1), coords(:,2), 'ko', 'MarkerFaceColor','k','MarkerSize',5);
            % Highlight
            if ~isempty(highlight_sites)
                plot(coords(highlight_sites,1), coords(highlight_sites,2), ...
                     'ro', 'MarkerFaceColor','r','MarkerSize',8);
            end
            axis equal; axis tight;
            xlabel('x'); ylabel('y');
        end

        % ---- Print geometry summary ---------------------------------
        function info(Nx, Ny, bcx, bcy)
            fprintf('\n  2D Square Lattice  %d × %d  =  %d sites\n', Nx, Ny, Nx*Ny);
            fprintf('  BC   : x=%s,  y=%s\n', upper(bcx), upper(bcy));
            fprintf('  Index: n(ix,iy) = (ix-1)*Ny + iy   (column-major)\n');
            fprintf('  Nambu block offsets:\n');
            Ns = Nx*Ny;
            fprintf('    c↑   : 1     .. %d\n',   Ns);
            fprintf('    c↓   : %d .. %d\n',   Ns+1,   2*Ns);
            fprintf('    c↑†  : %d .. %d\n',  2*Ns+1,  3*Ns);
            fprintf('    c↓†  : %d .. %d\n',  3*Ns+1,  4*Ns);
            fprintf('\n');
        end
    end
end
