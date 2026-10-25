classdef Lattice < UQMaT.Models.Base
    % UQMaT.Models.Lattice3D.Cubic.Lattice
    %   Simple cubic 3D tight-binding lattice.
    %
    %   Site index: n(ix,iy,iz) = (ix-1)*Ny*Nz + (iy-1)*Nz + iz
    %   Total sites: Ns = Nx * Ny * Nz
    %
    %   Lattice vectors: a1=[1,0,0], a2=[0,1,0], a3=[0,0,1]
    %   6 NN bonds per site: ±x, ±y, ±z
    %
    %   H sizes:
    %     normal  →  Ns × Ns
    %     spinful →  2Ns × 2Ns
    %     sc      →  4Ns × 4Ns
    %
    %   Anisotropic hopping: t = [tx, ty, tz]
    %
    %   Dispersion: E(k) = -2tx*cos(kx) - 2ty*cos(ky) - 2tz*cos(kz) - mu
    %
    %   Usage:
    %     p   = UQMaT.Params('system','cubic','Nx',6,'Ny',6,'Nz',6,'t',[1,1,0.5])
    %     sys = UQMaT.System(p)
    %     [kd,Ek,kl,kt] = UQMaT.Models.Lattice3D.Cubic.BC.bands(sys)

    %  Citation
    %  --------
    %  Simple cubic tight-binding: Ashcroft N.W., Mermin N.D., "Solid State
    %    Physics", Holt (1976). [Chapter 10]
    %  3D topological SC: Fu L., Berg E., PRL 105, 097001 (2010).
    %    [odd-parity SC on cubic; Cu_x Bi_2 Se_3]
    methods
        function obj = Lattice(p)
            obj.params = p;
            obj.Ndim   = p.Nx * p.Ny * p.Nz;
            obj.H      = obj.buildH();
        end

        function s = description(obj)
            p  = obj.params;
            Ns = p.Nx*p.Ny*p.Nz;
            s  = sprintf('Cubic3D %d×%d×%d=%d sites  bc=%s', ...
                p.Nx, p.Ny, p.Nz, Ns, upper(p.bc));
        end

        function H = buildH(obj)
            p   = obj.params;
            Nx  = p.Nx; Ny = p.Ny; Nz = p.Nz;
            Ns  = Nx*Ny*Nz;
            bdg = UQMaT.SC.BdG;
            pbc = strcmp(p.bc,'pbc');

            % Disorder and OPS
            [D_s_ops, D_t_ops, h_AM, h_CDW] = UQMaT.OrderParameterSymmetry.OPS.build(p);
            V_dis = diag(UQMaT.Disorder.Disorder.onsite(p, Ns));

            % Build kinetic h0
            h0 = obj.hopping_matrix(Nx,Ny,Nz,p.tx,p.ty,p.tz,p.t2,pbc) ...
                 - p.mu*eye(Ns) + p.onsite_val*eye(Ns) + V_dis + h_CDW;

            if p.sc
                H_N = kron(eye(2), h0);
                if any(p.B ~= 0), H_N = H_N + bdg.zeeman(p.B,Ns); end
                H_N = H_N + h_AM;
                D   = D_s_ops + D_t_ops;
                H   = bdg.assemble(H_N, D);
            elseif p.spin
                s   = bdg.sigma();
                H   = kron(eye(2), h0);
                if any(p.B ~= 0), H = H + bdg.zeeman(p.B,Ns); end
                H   = H + h_AM;
            else
                H   = h0;
            end
        end
    end

    methods (Static)

        function n = site_index(ix,iy,iz,Ny,Nz)
            n = (ix-1)*Ny*Nz + (iy-1)*Nz + iz;
        end

        function T = hopping_matrix(Nx,Ny,Nz,tx,ty,tz,t2,pbc)
            Ns = Nx*Ny*Nz;
            T  = zeros(Ns);
            id = @(x,y,z) UQMaT.Models.Lattice3D.Cubic.Lattice.idx_wrap(x,y,z,Nx,Ny,Nz,pbc);

            for ix=1:Nx; for iy=1:Ny; for iz=1:Nz
                n0 = (ix-1)*Ny*Nz + (iy-1)*Nz + iz;
                % x-bonds
                m=id(ix+1,iy,iz); if ~isnan(m), T(n0,m)=T(n0,m)-tx; T(m,n0)=T(m,n0)-tx; end
                % y-bonds
                m=id(ix,iy+1,iz); if ~isnan(m), T(n0,m)=T(n0,m)-ty; T(m,n0)=T(m,n0)-ty; end
                % z-bonds
                m=id(ix,iy,iz+1); if ~isnan(m), T(n0,m)=T(n0,m)-tz; T(m,n0)=T(m,n0)-tz; end
                % NNN (if t2 nonzero): face diagonals
                if t2 ~= 0
                    for dnbr = {[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]}
                        d=dnbr{1};
                        m=id(ix+d(1),iy+d(2),iz+d(3));
                        if ~isnan(m), T(n0,m)=T(n0,m)-t2; T(m,n0)=T(m,n0)-t2; end
                    end
                end
            end; end; end
        end

    end

    methods (Static, Access=private)
        function n = idx_wrap(ix,iy,iz,Nx,Ny,Nz,pbc)
            if pbc
                ix=mod(ix-1,Nx)+1; iy=mod(iy-1,Ny)+1; iz=mod(iz-1,Nz)+1;
            end
            if ix<1||ix>Nx||iy<1||iy>Ny||iz<1||iz>Nz, n=NaN; return; end
            n = (ix-1)*Ny*Nz + (iy-1)*Nz + iz;
        end
    end
end
