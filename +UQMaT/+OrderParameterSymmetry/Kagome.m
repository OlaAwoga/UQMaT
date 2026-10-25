classdef Kagome
    % UQMaT.OrderParameterSymmetry.Kagome
    %   Order parameter bond structures for the 2D kagome lattice.
    %
    %   3 sublattices (A,B,C) per unit cell; Ns = 3*Nx*Ny sites.
    %   NN bonds per unit cell: A-B, A-C, B-C (intra) + 3 inter-cell.
    %   NNN bonds: same-sublattice, 4 per site.
    %
    %   Symmetries:
    %     onsite      : Δ(i) = 1 (s-wave, uniform)
    %     extended_s  : equal weight on all 6 NN bonds per cell
    %     px          : antisymmetric along x-direction bonds
    %     py          : antisymmetric along y-direction bonds
    %     chiral_p    : [1, e^{2πi/3}, ...] on NN bonds (breaks TRS)
    %     dx2y2       : sign-alternating on NN bonds
    %     d_chiral    : NNN bonds with imaginary phase (d+id)
    %
    %   Usage:
    %     sym = UQMaT.OrderParameterSymmetry.Kagome.get('chiral_p', Nx, Ny, bc)

    %  Citation
    %  --------
    %  Kagome SC: Ko W.-H. et al., PRB 79, 214502 (2009).
    %  Chiral SC on kagome: Mazin I.I. et al., Nat. Comm. 5, 4261 (2014).
    %  Sigrist M., Ueda K., Rev. Mod. Phys. 63, 239 (1991).
    methods (Static)

        function sym = get(label, Nx, Ny, bc)
            if nargin < 4, bc = 'obc'; end
            label = lower(strtrim(label));
            switch label
                case 'onsite'
                    sym = UQMaT.OrderParameterSymmetry.Kagome.onsite_(Nx,Ny);
                case 'extended_s'
                    sym = UQMaT.OrderParameterSymmetry.Kagome.nn_matrix_(Nx,Ny,bc, ...
                          ones(6,1)/sqrt(6));
                case 'px'
                    % Antisymmetric on x-direction bonds (B-A inter-cell)
                    w = 1/sqrt(2);
                    wts = [0; 0; 0; +w; -w; 0];  % only x-direction bonds
                    sym = UQMaT.OrderParameterSymmetry.Kagome.nn_matrix_(Nx,Ny,bc,wts);
                case 'py'
                    w = 1/sqrt(2);
                    wts = [0; 0; 0; 0; +w; -w];  % y-direction bonds
                    sym = UQMaT.OrderParameterSymmetry.Kagome.nn_matrix_(Nx,Ny,bc,wts);
                case 'chiral_p'
                    % [1, e^{2πi/3}, e^{4πi/3}, ...]/√6 on NN bonds
                    ph  = exp(1i*2*pi/3);
                    wts = [1; ph; ph^2; ph^3; ph^4; ph^5] / sqrt(6);
                    sym = UQMaT.OrderParameterSymmetry.Kagome.nn_matrix_(Nx,Ny,bc,wts);
                case 'dx2y2'
                    % Alternating signs on NN bonds
                    wts = [+1;-1;+1;-1;+1;-1] / sqrt(6);
                    sym = UQMaT.OrderParameterSymmetry.Kagome.nn_matrix_(Nx,Ny,bc,wts);
                case 'd_chiral'
                    % NNN bonds with imaginary phase (d+id)
                    sym = UQMaT.OrderParameterSymmetry.Kagome.nnn_matrix_(Nx,Ny,bc, ...
                          1i*ones(4,1)/2);
                otherwise
                    error('UQMaT:OPS:Kagome:unknown', ...
                        ['Unknown symmetry "%s".\n' ...
                         'Available: onsite, extended_s, px, py, chiral_p, dx2y2, d_chiral'], ...
                        label);
            end
            sym.label  = label;
            sym.system = 'kagome';
        end

        function list = available()
            list = {'onsite','extended_s','px','py','chiral_p','dx2y2','d_chiral'};
        end

    end

    methods (Static, Access = private)

        function sym = onsite_(Nx, Ny)
            Ns = 3*Nx*Ny;
            sym.D            = eye(Ns);
            sym.site_weights = ones(Ns,1);
            % Sublattice: A=+1, B=-1, C=+1 (kagome has 3-site cell)
            sub = zeros(Ns,1);
            for ix=1:Nx; for iy=1:Ny
                n = @(s) UQMaT.Models.Lattice2D.Kagome.Geometry.idx(ix,iy,s,Ny);
                sub(n(1))=+1; sub(n(2))=-1; sub(n(3))=+1;
            end; end
            sym.sublattice = sub;
        end

        function n = kidx(ix,iy,sub,Nx,Ny,pbc)
            if pbc, ix=mod(ix-1,Nx)+1; iy=mod(iy-1,Ny)+1; end
            if ix<1||ix>Nx||iy<1||iy>Ny, n=NaN; return; end
            n = UQMaT.Models.Lattice2D.Kagome.Geometry.idx(ix,iy,sub,Ny);
        end

        function sym = nn_matrix_(Nx, Ny, bc, wts)
            % Build Ns×Ns bond matrix with weights on 6 NN bonds per cell.
            % Bond order: AB-intra, AC-intra, BC-intra, BA-inter(x),
            %             CA-inter(y), CB-inter(diag)
            Ns  = 3*Nx*Ny;
            D   = zeros(Ns);
            pbc = strcmp(bc,'pbc');
            id  = @(ix,iy,s) UQMaT.OrderParameterSymmetry.Kagome.kidx(ix,iy,s,Nx,Ny,pbc);

            for ix=1:Nx; for iy=1:Ny
                nA=id(ix,iy,1); nB=id(ix,iy,2); nC=id(ix,iy,3);
                % Intra-cell
                if ~isnan(nA)&&~isnan(nB), D(nA,nB)=D(nA,nB)+wts(1); D(nB,nA)=D(nB,nA)+conj(wts(1)); end
                if ~isnan(nA)&&~isnan(nC), D(nA,nC)=D(nA,nC)+wts(2); D(nC,nA)=D(nC,nA)+conj(wts(2)); end
                if ~isnan(nB)&&~isnan(nC), D(nB,nC)=D(nB,nC)+wts(3); D(nC,nB)=D(nC,nB)+conj(wts(3)); end
                % Inter-cell
                m=id(ix+1,iy,1);   if ~isnan(nB)&&~isnan(m), D(nB,m)=D(nB,m)+wts(4); D(m,nB)=D(m,nB)+conj(wts(4)); end
                m=id(ix,iy+1,1);   if ~isnan(nC)&&~isnan(m), D(nC,m)=D(nC,m)+wts(5); D(m,nC)=D(m,nC)+conj(wts(5)); end
                m=id(ix-1,iy+1,2); if ~isnan(nC)&&~isnan(m), D(nC,m)=D(nC,m)+wts(6); D(m,nC)=D(m,nC)+conj(wts(6)); end
            end; end
            sym.D = D;
            sym.site_weights = zeros(Ns,1);
            sym.sublattice   = UQMaT.OrderParameterSymmetry.Kagome.onsite_(Nx,Ny).sublattice;
        end

        function sym = nnn_matrix_(Nx, Ny, bc, wts)
            % NNN same-sublattice bonds with weights wts (4×1)
            Ns  = 3*Nx*Ny;
            D   = zeros(Ns);
            pbc = strcmp(bc,'pbc');
            id  = @(ix,iy,s) UQMaT.OrderParameterSymmetry.Kagome.kidx(ix,iy,s,Nx,Ny,pbc);
            dirs= [+1,0;-1,0;0,+1;0,-1];
            for ix=1:Nx; for iy=1:Ny
                for sub=1:3
                    n0=id(ix,iy,sub);
                    if isnan(n0), continue; end
                    for k=1:4
                        n1=id(ix+dirs(k,1),iy+dirs(k,2),sub);
                        if ~isnan(n1), D(n0,n1)=D(n0,n1)+wts(k); end
                    end
                end
            end; end
            sym.D = D;
            sym.site_weights = zeros(Ns,1);
            sym.sublattice   = UQMaT.OrderParameterSymmetry.Kagome.onsite_(Nx,Ny).sublattice;
        end

    end
end
