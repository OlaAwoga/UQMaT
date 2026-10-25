classdef Square2D
    % UQMaT.OrderParameterSymmetry.Square2D
    %   Order parameter bond structures for the 2D square lattice.
    %
    %   Site index: n(ix,iy) = (ix-1)*Ny + iy  (column-major)
    %
    %   NN bonds (4 per site):  (x-1,y), (x+1,y), (x,y-1), (x,y+1)
    %   NNN bonds (4 per site): (x+1,y+1), (x-1,y+1), (x+1,y-1), (x-1,y-1)
    %
    %   Symmetries:
    %     onsite      : Δ(i) = 1  (s-wave, uniform)
    %     extended_s  : [+1,+1,+1,+1]/2  on 4 NN bonds  (cos kx + cos ky)
    %     px          : [+1,-1, 0, 0]/√2 on NN (x-bonds) (sin kx)
    %     py          : [ 0, 0,+1,-1]/√2 on NN (y-bonds) (sin ky)
    %     chiral_p    : [+1,-1,+i,-i]/2  on 4 NN bonds   (sin kx + i sin ky)
    %     dx2y2       : [+1,+1,-1,-1]/2  on 4 NN bonds   (cos kx - cos ky)
    %     dxy         : [+1,-1,-1,+1]/2  on 4 NNN bonds  (sin kx sin ky)
    %
    %   Usage:
    %     sym = UQMaT.OrderParameterSymmetry.Square2D.get('dx2y2', Nx, Ny, bcx, bcy)
    %     D   = Delta * sym.D   % Ns×Ns bond matrix

    %  Citation
    %  --------
    %  d-wave SC: Zhang F.C., Rice T.M., PRB 37, 3759 (1988).
    %  Chiral p-wave: Read N., Green D., PRB 61, 10267 (2000).
    %  Sigrist M., Ueda K., Rev. Mod. Phys. 63, 239 (1991). [gap symmetry]
    methods (Static)

        function sym = get(label, Nx, Ny, bcx, bcy)
            if nargin < 4, bcx = 'obc'; end
            if nargin < 5, bcy = 'obc'; end
            label = lower(strtrim(label));
            switch label
                case 'onsite'
                    sym = UQMaT.OrderParameterSymmetry.Square2D.onsite_(Nx,Ny);
                case 'extended_s'
                    sym = UQMaT.OrderParameterSymmetry.Square2D.extended_s_(Nx,Ny,bcx,bcy);
                case 'px'
                    sym = UQMaT.OrderParameterSymmetry.Square2D.px_(Nx,Ny,bcx,bcy);
                case 'py'
                    sym = UQMaT.OrderParameterSymmetry.Square2D.py_(Nx,Ny,bcx,bcy);
                case 'chiral_p'
                    sym = UQMaT.OrderParameterSymmetry.Square2D.chiral_p_(Nx,Ny,bcx,bcy);
                case 'dx2y2'
                    sym = UQMaT.OrderParameterSymmetry.Square2D.dx2y2_(Nx,Ny,bcx,bcy);
                case 'dxy'
                    sym = UQMaT.OrderParameterSymmetry.Square2D.dxy_(Nx,Ny,bcx,bcy);
                otherwise
                    error('UQMaT:OPS:Square2D:unknown', ...
                        ['Unknown symmetry "%s".\n' ...
                         'Available: onsite, extended_s, px, py, ' ...
                         'chiral_p, dx2y2, dxy'], label);
            end
            sym.label  = label;
            sym.system = 'square';
        end

        function list = available()
            list = {'onsite','extended_s','px','py','chiral_p','dx2y2','dxy'};
        end

    end

    methods (Static, Access = private)

        function n = idx(ix, iy, Nx, Ny)
            % Column-major site index: n(ix,iy) = (ix-1)*Ny + iy
            % Returns NaN for out-of-bounds (OBC)
            if ix < 1 || ix > Nx || iy < 1 || iy > Ny
                n = NaN;
            else
                n = (ix-1)*Ny + iy;
            end
        end

        function n = idx_pbc(ix, iy, Nx, Ny, bcx, bcy)
            % Site index with periodic wrapping
            if strcmp(bcx,'pbc'), ix = mod(ix-1,Nx)+1; end
            if strcmp(bcy,'pbc'), iy = mod(iy-1,Ny)+1; end
            n = UQMaT.OrderParameterSymmetry.Square2D.idx(ix,iy,Nx,Ny);
        end

        function sym = onsite_(Nx, Ny)
            Ns = Nx*Ny;
            sym.D            = eye(Ns);
            sym.site_weights = ones(Ns,1);
            % Checkerboard sublattice for AM
            [IX,IY] = meshgrid(1:Nx,1:Ny);
            sub     = reshape((-1).^(IX+IY)', Ns, 1);
            sym.sublattice = sub;
        end

        function sym = nn_bond_matrix_(Nx, Ny, bcx, bcy, wx, wy)
            % Build Ns×Ns bond matrix for NN bonds with weights:
            %   wx(1) = weight for (ix-1,iy) bond
            %   wx(2) = weight for (ix+1,iy) bond
            %   wy(1) = weight for (ix,iy-1) bond
            %   wy(2) = weight for (ix,iy+1) bond
            Ns = Nx*Ny;
            D  = zeros(Ns);
            id = @(ix,iy) UQMaT.OrderParameterSymmetry.Square2D.idx_pbc(ix,iy,Nx,Ny,bcx,bcy);

            for ix = 1:Nx
                for iy = 1:Ny
                    n = (ix-1)*Ny + iy;
                    % x-direction bonds
                    m = id(ix-1,iy); if ~isnan(m), D(n,m) = wx(1); end
                    m = id(ix+1,iy); if ~isnan(m), D(n,m) = wx(2); end
                    % y-direction bonds
                    m = id(ix,iy-1); if ~isnan(m), D(n,m) = wy(1); end
                    m = id(ix,iy+1); if ~isnan(m), D(n,m) = wy(2); end
                end
            end
            sym.D = D;
            sym.site_weights = zeros(Ns,1);
            [IX,IY] = meshgrid(1:Nx,1:Ny);
            sym.sublattice = reshape((-1).^(IX+IY)', Ns, 1);
        end

        function sym = nnn_bond_matrix_(Nx, Ny, bcx, bcy, weights)
            % NNN bonds: (±1,±1) with weights [w++, w-+, w+-, w--]
            %   w++ : (ix+1,iy+1)
            %   w-+ : (ix-1,iy+1)
            %   w+- : (ix+1,iy-1)
            %   w-- : (ix-1,iy-1)
            Ns = Nx*Ny;
            D  = zeros(Ns);
            id = @(ix,iy) UQMaT.OrderParameterSymmetry.Square2D.idx_pbc(ix,iy,Nx,Ny,bcx,bcy);
            dirs = [+1,+1; -1,+1; +1,-1; -1,-1];
            for ix = 1:Nx
                for iy = 1:Ny
                    n = (ix-1)*Ny + iy;
                    for k = 1:4
                        m = id(ix+dirs(k,1), iy+dirs(k,2));
                        if ~isnan(m), D(n,m) = weights(k); end
                    end
                end
            end
            sym.D = D;
            sym.site_weights = zeros(Ns,1);
            [IX,IY] = meshgrid(1:Nx,1:Ny);
            sym.sublattice = reshape((-1).^(IX+IY)', Ns, 1);
        end

        function sym = extended_s_(Nx, Ny, bcx, bcy)
            % Extended s-wave: [+1,+1,+1,+1]/2 on 4 NN bonds
            % Fourier: cos(kx) + cos(ky)
            w = 1/2;
            sym = UQMaT.OrderParameterSymmetry.Square2D.nn_bond_matrix_(...
                Nx,Ny,bcx,bcy, [w,w], [w,w]);
        end

        function sym = px_(Nx, Ny, bcx, bcy)
            % p_x: antisymmetric in x, zero in y
            % weights: [+1,-1,0,0]/√2
            w = 1/sqrt(2);
            sym = UQMaT.OrderParameterSymmetry.Square2D.nn_bond_matrix_(...
                Nx,Ny,bcx,bcy, [+w,-w], [0,0]);
        end

        function sym = py_(Nx, Ny, bcx, bcy)
            % p_y: antisymmetric in y, zero in x
            % weights: [0,0,+1,-1]/√2
            w = 1/sqrt(2);
            sym = UQMaT.OrderParameterSymmetry.Square2D.nn_bond_matrix_(...
                Nx,Ny,bcx,bcy, [0,0], [+w,-w]);
        end

        function sym = chiral_p_(Nx, Ny, bcx, bcy)
            % Chiral p-wave: p_x + i*p_y
            % weights: [+1,-1,+i,-i]/2
            % Breaks TRS, non-zero Chern number (p+ip SC)
            w = 1/2;
            sym = UQMaT.OrderParameterSymmetry.Square2D.nn_bond_matrix_(...
                Nx,Ny,bcx,bcy, [+w,-w], [+1i*w,-1i*w]);
        end

        function sym = dx2y2_(Nx, Ny, bcx, bcy)
            % d_{x²-y²}: [+1,+1,-1,-1]/2 on NN bonds
            % Fourier: cos(kx) - cos(ky)
            w = 1/2;
            sym = UQMaT.OrderParameterSymmetry.Square2D.nn_bond_matrix_(...
                Nx,Ny,bcx,bcy, [+w,+w], [-w,-w]);
        end

        function sym = dxy_(Nx, Ny, bcx, bcy)
            % d_{xy}: NNN bonds [+1,-1,-1,+1]/2
            % Fourier: sin(kx)sin(ky)
            w = 1/2;
            sym = UQMaT.OrderParameterSymmetry.Square2D.nnn_bond_matrix_(...
                Nx,Ny,bcx,bcy, [+w,-w,-w,+w]);
        end

    end
end
