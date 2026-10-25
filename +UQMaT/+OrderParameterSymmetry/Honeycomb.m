classdef Honeycomb
    % UQMaT.OrderParameterSymmetry.Honeycomb
    %   Order parameter bond structures for the honeycomb lattice.
    %
    %   Honeycomb has two sublattices A and B per unit cell.
    %   Site index: n(ix,iy,sub) where sub ∈ {1=A, 2=B}
    %   Total sites: Ns = 2*Nx*Ny
    %
    %   NN bonds (3 per site): δ₁=(0,0), δ₂=(1,0), δ₃=(0,1) from A to B
    %   NNN bonds (6 per site): connect same-sublattice sites
    %     A-A: (±1,0), (0,±1), (1,-1), (-1,+1)
    %     B-B: same
    %
    %   Symmetries:
    %     onsite      : Δ(i) = 1 (s-wave, uniform)
    %     extended_s  : [+1,+1,+1]/√3 on 3 NN bonds
    %     px          : Re[e^{ik·δ}] weighting on NN bonds
    %     py          : Im[e^{ik·δ}] weighting on NN bonds
    %     chiral_p    : [1, e^{2πi/3}, e^{4πi/3}]/√3  (breaks TRS)
    %     dx2y2       : NNN bonds, d-wave pattern
    %     d_chiral    : NNN bonds, d+id (Chern insulator SC)
    %
    %   Usage:
    %     sym = UQMaT.OrderParameterSymmetry.Honeycomb.get('chiral_p', Nx, Ny, bc)

    %  Citation
    %  --------
    %  Honeycomb lattice: Wallace P.R., Phys. Rev. 71, 622 (1947). [graphene TB]
    %  Chiral SC on honeycomb: Black-Schaffer A.M., Doniach S., PRB 75, 134512 (2007).
    %  d+id SC: Nandkishore R. et al., Nat. Phys. 8, 158 (2012).
    methods (Static)

        function sym = get(label, Nx, Ny, bc)
            if nargin < 4, bc = 'obc'; end
            label = lower(strtrim(label));
            switch label
                case 'onsite'
                    sym = UQMaT.OrderParameterSymmetry.Honeycomb.onsite_(Nx,Ny);
                case 'extended_s'
                    sym = UQMaT.OrderParameterSymmetry.Honeycomb.extended_s_(Nx,Ny,bc);
                case 'px'
                    sym = UQMaT.OrderParameterSymmetry.Honeycomb.px_(Nx,Ny,bc);
                case 'py'
                    sym = UQMaT.OrderParameterSymmetry.Honeycomb.py_(Nx,Ny,bc);
                case 'chiral_p'
                    sym = UQMaT.OrderParameterSymmetry.Honeycomb.chiral_p_(Nx,Ny,bc);
                case 'dx2y2'
                    sym = UQMaT.OrderParameterSymmetry.Honeycomb.dx2y2_(Nx,Ny,bc);
                case 'd_chiral'
                    sym = UQMaT.OrderParameterSymmetry.Honeycomb.d_chiral_(Nx,Ny,bc);
                otherwise
                    error('UQMaT:OPS:Honeycomb:unknown', ...
                        ['Unknown symmetry "%s".\n' ...
                         'Available: onsite, extended_s, px, py, ' ...
                         'chiral_p, dx2y2, d_chiral'], label);
            end
            sym.label  = label;
            sym.system = 'honeycomb';
        end

        function list = available()
            list = {'onsite','extended_s','px','py','chiral_p','dx2y2','d_chiral'};
        end

        function n = site_index(ix, iy, sub, Nx, Ny)
            % n(ix,iy,sub) = 2*((ix-1)*Ny + (iy-1)) + sub
            % sub: 1=A, 2=B
            n = 2*((ix-1)*Ny + (iy-1)) + sub;
        end

    end

    methods (Static, Access = private)

        function n = idx(ix, iy, sub, Nx, Ny, bc)
            % Honeycomb site index with periodic wrapping
            if strcmp(bc,'pbc')
                ix = mod(ix-1,Nx)+1;
                iy = mod(iy-1,Ny)+1;
            end
            if ix < 1 || ix > Nx || iy < 1 || iy > Ny
                n = NaN;
            else
                n = UQMaT.OrderParameterSymmetry.Honeycomb.site_index(ix,iy,sub,Nx,Ny);
            end
        end

        function sym = onsite_(Nx, Ny)
            Ns = 2*Nx*Ny;
            sym.D            = eye(Ns);
            sym.site_weights = ones(Ns,1);
            % Sublattice: A sites = +1, B sites = -1
            sub = zeros(Ns,1);
            for ix=1:Nx; for iy=1:Ny
                sub(UQMaT.OrderParameterSymmetry.Honeycomb.site_index(ix,iy,1,Nx,Ny)) = +1;
                sub(UQMaT.OrderParameterSymmetry.Honeycomb.site_index(ix,iy,2,Nx,Ny)) = -1;
            end; end
            sym.sublattice = sub;
        end

        function sym = nn_bond_matrix_(Nx, Ny, bc, weights)
            % NN bonds from A→B:
            %   δ₁: same cell       (ix,iy)   A→B  weight(1)
            %   δ₂: (ix+1,iy)       A→B       weight(2)
            %   δ₃: (ix,iy+1)       A→B       weight(3)
            % Matrix is NOT symmetric (A→B ≠ B→A for chiral)
            Ns = 2*Nx*Ny;
            D  = zeros(Ns);
            id = @(ix,iy,s) UQMaT.OrderParameterSymmetry.Honeycomb.idx(ix,iy,s,Nx,Ny,bc);

            for ix = 1:Nx
                for iy = 1:Ny
                    nA = id(ix,iy,1);
                    % δ₁: (ix,iy) A→B
                    nB = id(ix,iy,2);
                    if ~isnan(nA) && ~isnan(nB)
                        D(nA,nB) = D(nA,nB) + weights(1);
                        D(nB,nA) = D(nB,nA) + conj(weights(1));
                    end
                    % δ₂: (ix+1,iy) A→B
                    nB = id(ix+1,iy,2);
                    if ~isnan(nA) && ~isnan(nB)
                        D(nA,nB) = D(nA,nB) + weights(2);
                        D(nB,nA) = D(nB,nA) + conj(weights(2));
                    end
                    % δ₃: (ix,iy+1) A→B
                    nB = id(ix,iy+1,2);
                    if ~isnan(nA) && ~isnan(nB)
                        D(nA,nB) = D(nA,nB) + weights(3);
                        D(nB,nA) = D(nB,nA) + conj(weights(3));
                    end
                end
            end
            sym.D            = D;
            sym.site_weights = zeros(Ns,1);
            % Sublattice pattern for AM
            sub = zeros(Ns,1);
            for ix=1:Nx; for iy=1:Ny
                sub(UQMaT.OrderParameterSymmetry.Honeycomb.site_index(ix,iy,1,Nx,Ny)) = +1;
                sub(UQMaT.OrderParameterSymmetry.Honeycomb.site_index(ix,iy,2,Nx,Ny)) = -1;
            end; end
            sym.sublattice = sub;
        end

        function sym = nnn_bond_matrix_(Nx, Ny, bc, wA, wB)
            % NNN bonds: same-sublattice hoppings
            % A-A NNN: (±1,0), (0,±1) relative
            % wA = 6-element weight vector for A sites (±1,0),(0,±1),(1,-1),(-1,1)
            % wB = same for B sites
            Ns   = 2*Nx*Ny;
            D    = zeros(Ns);
            id   = @(ix,iy,s) UQMaT.OrderParameterSymmetry.Honeycomb.idx(ix,iy,s,Nx,Ny,bc);
            % NNN directions for A sublattice
            dirs = [+1,0; -1,0; 0,+1; 0,-1; +1,-1; -1,+1];
            for ix = 1:Nx
                for iy = 1:Ny
                    nA0 = id(ix,iy,1);
                    nB0 = id(ix,iy,2);
                    for k = 1:6
                        dx = dirs(k,1); dy = dirs(k,2);
                        % A-A
                        nA1 = id(ix+dx,iy+dy,1);
                        if ~isnan(nA0) && ~isnan(nA1)
                            D(nA0,nA1) = D(nA0,nA1) + wA(k);
                        end
                        % B-B
                        nB1 = id(ix+dx,iy+dy,2);
                        if ~isnan(nB0) && ~isnan(nB1)
                            D(nB0,nB1) = D(nB0,nB1) + wB(k);
                        end
                    end
                end
            end
            sym.D            = D;
            sym.site_weights = zeros(Ns,1);
            sub = zeros(Ns,1);
            for ix=1:Nx; for iy=1:Ny
                sub(UQMaT.OrderParameterSymmetry.Honeycomb.site_index(ix,iy,1,Nx,Ny)) = +1;
                sub(UQMaT.OrderParameterSymmetry.Honeycomb.site_index(ix,iy,2,Nx,Ny)) = -1;
            end; end
            sym.sublattice = sub;
        end

        function sym = extended_s_(Nx, Ny, bc)
            % Extended s-wave: equal weight [1,1,1]/√3 on all 3 NN bonds
            w = 1/sqrt(3);
            sym = UQMaT.OrderParameterSymmetry.Honeycomb.nn_bond_matrix_(Nx,Ny,bc,[w,w,w]);
        end

        function sym = px_(Nx, Ny, bc)
            % p_x: real part of e^{ik·δ} weighting
            % weights ∝ Re[1, e^{i2π/3}, e^{i4π/3}] = [1, -1/2, -1/2]/√(3/2)
            w = sqrt(2/3);
            wts = w * [1, -0.5, -0.5];
            sym = UQMaT.OrderParameterSymmetry.Honeycomb.nn_bond_matrix_(Nx,Ny,bc,wts);
        end

        function sym = py_(Nx, Ny, bc)
            % p_y: imaginary part of e^{ik·δ} weighting
            % weights ∝ Im[1, e^{i2π/3}, e^{i4π/3}] = [0, √3/2, -√3/2]
            w = 1/sqrt(2);
            wts = w * [0, sqrt(3)/2, -sqrt(3)/2] * (2/sqrt(3));
            sym = UQMaT.OrderParameterSymmetry.Honeycomb.nn_bond_matrix_(Nx,Ny,bc,wts);
        end

        function sym = chiral_p_(Nx, Ny, bc)
            % Chiral p-wave (p+ip): [1, e^{2πi/3}, e^{4πi/3}]/√3
            % Breaks TRS → non-zero Chern number
            w = 1/sqrt(3);
            wts = w * [1, exp(1i*2*pi/3), exp(1i*4*pi/3)];
            sym = UQMaT.OrderParameterSymmetry.Honeycomb.nn_bond_matrix_(Nx,Ny,bc,wts);
        end

        function sym = dx2y2_(Nx, Ny, bc)
            % d_{x²-y²}: NNN bonds with alternating signs
            % Haldane-like pattern: wA = +1 for CCW, wB = -1
            w  = 1/sqrt(6);
            wA = w * [+1,-1,+1,-1,+1,-1];
            wB = w * [-1,+1,-1,+1,-1,+1];
            sym = UQMaT.OrderParameterSymmetry.Honeycomb.nnn_bond_matrix_(Nx,Ny,bc,wA,wB);
        end

        function sym = d_chiral_(Nx, Ny, bc)
            % d+id: NNN bonds with complex Haldane-like pattern
            % Uniform imaginary phase on NNN → topological SC
            w  = 1/sqrt(6);
            ph = exp(1i*pi/2);   % i
            wA = w * ph * ones(1,6);
            wB = w * conj(ph) * ones(1,6);
            sym = UQMaT.OrderParameterSymmetry.Honeycomb.nnn_bond_matrix_(Nx,Ny,bc,wA,wB);
        end

    end
end
