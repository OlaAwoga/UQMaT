classdef Lattice < UQMaT.Models.Base
    % UQMaT.Models.Lattice3D.Diamond.Lattice
    %   Diamond cubic lattice: FCC Bravais lattice with 2-site basis.
    %
    %   Basis: A at (0,0,0), B at (a/4)(1,1,1)
    %   FCC lattice vectors:
    %     a1 = a/2 * [0,1,1]
    %     a2 = a/2 * [1,0,1]
    %     a3 = a/2 * [1,1,0]
    %
    %   Each A site has 4 NN B sites; each B site has 4 NN A sites.
    %   Total sites: Ns = 2 * Nx * Ny * Nz
    %   Site index: n(ix,iy,iz,sub) = 2*((ix-1)*Ny*Nz + (iy-1)*Nz + (iz-1)) + sub
    %               sub: 1=A, 2=B
    %
    %   Spectrum: 2 bands, gapped for any nonzero t.
    %             At Γ: E = ±4t (bonding/antibonding)
    %
    %   NN vectors from A to B (in units of a/4):
    %     d1 = (+1,+1,+1)    d2 = (+1,-1,-1)
    %     d3 = (-1,+1,-1)    d4 = (-1,-1,+1)
    %
    %   Usage:
    %     p   = UQMaT.Params('system','diamond','Nx',4,'Ny',4,'Nz',4,'t',1)
    %     sys = UQMaT.System(p)

    %  Citation
    %  --------
    %  Diamond structure: Ashcroft N.W., Mermin N.D., "Solid State Physics",
    %    Holt (1976). [Chapter 4: FCC basis, diamond bond geometry]
    %  Harrison W.A., "Electronic Structure and the Properties of Solids",
    %    Freeman (1980). [Chapter 4: sp3 hybridisation]
    methods
        function obj = Lattice(p)
            obj.params = p;
            obj.Ndim   = 2 * p.Nx * p.Ny * p.Nz;
            obj.H      = obj.buildH();
        end

        function s = description(obj)
            p  = obj.params;
            Ns = 2*p.Nx*p.Ny*p.Nz;
            s  = sprintf('Diamond3D %d×%d×%d=%d sites  bc=%s', ...
                p.Nx, p.Ny, p.Nz, Ns, upper(p.bc));
        end

        function H = buildH(obj)
            p   = obj.params;
            Nx  = p.Nx; Ny = p.Ny; Nz = p.Nz;
            Ns  = 2*Nx*Ny*Nz;
            bdg = UQMaT.SC.BdG;
            pbc = strcmp(p.bc,'pbc');

            [D_s_ops, D_t_ops, h_AM, h_CDW] = UQMaT.OrderParameterSymmetry.OPS.build(p);
            V_dis = diag(UQMaT.Disorder.Disorder.onsite(p, Ns));

            h0 = obj.hopping_matrix(Nx,Ny,Nz,p.t,pbc) ...
                 - p.mu*eye(Ns) + p.onsite_val*eye(Ns) + V_dis + h_CDW;

            if p.sc
                H_N = kron(eye(2), h0);
                if any(p.B~=0), H_N = H_N + bdg.zeeman(p.B,Ns); end
                H_N = H_N + h_AM;
                H   = bdg.assemble(H_N, D_s_ops + D_t_ops);
            elseif p.spin
                H = kron(eye(2), h0);
                if any(p.B~=0), H = H + bdg.zeeman(p.B,Ns); end
                H = H + h_AM;
            else
                H = h0;
            end
        end
    end

    methods (Static)

        function n = site_index(ix,iy,iz,sub,Ny,Nz)
            % n = 2*((ix-1)*Ny*Nz + (iy-1)*Nz + (iz-1)) + sub
            n = 2*((ix-1)*Ny*Nz + (iy-1)*Nz + (iz-1)) + sub;
        end

        function T = hopping_matrix(Nx,Ny,Nz,t,pbc)
            % NN vectors from A to neighbouring B cells (in unit cell coords):
            % d1=(0,0,0), d2=(1,0,0), d3=(0,1,0), d4=(0,0,1)
            % (A in cell (ix,iy,iz) connects to B in cell (ix,iy,iz),
            %  (ix+1,iy,iz), (ix,iy+1,iz), (ix,iy,iz+1))
            Ns  = 2*Nx*Ny*Nz;
            T   = zeros(Ns);
            id  = @(x,y,z,s) UQMaT.Models.Lattice3D.Diamond.Lattice.idx_wrap(x,y,z,s,Nx,Ny,Nz,pbc);
            NN_B_offsets = [0,0,0; 1,0,0; 0,1,0; 0,0,1];

            for ix=1:Nx; for iy=1:Ny; for iz=1:Nz
                nA = id(ix,iy,iz,1);
                for k=1:4
                    dx=NN_B_offsets(k,1); dy=NN_B_offsets(k,2); dz=NN_B_offsets(k,3);
                    nB = id(ix+dx,iy+dy,iz+dz,2);
                    if ~isnan(nA)&&~isnan(nB)
                        T(nA,nB)=T(nA,nB)-t; T(nB,nA)=T(nB,nA)-t;
                    end
                end
            end; end; end
        end

    end

    methods (Static, Access=private)
        function n = idx_wrap(ix,iy,iz,sub,Nx,Ny,Nz,pbc)
            if pbc
                ix=mod(ix-1,Nx)+1; iy=mod(iy-1,Ny)+1; iz=mod(iz-1,Nz)+1;
            end
            if ix<1||ix>Nx||iy<1||iy>Ny||iz<1||iz>Nz, n=NaN; return; end
            n = UQMaT.Models.Lattice3D.Diamond.Lattice.site_index(ix,iy,iz,sub,Ny,Nz);
        end
    end
end
