classdef Geometry
    % UQMaT.Models.Lattice2D.Square.Geometry
    %   Site indexing and k-path utilities for the 2D square lattice.
    %
    %   Site index: n(ix,iy) = (ix-1)*Ny + iy  (column-major)
    %   BZ high-symmetry path: Γ→X→M→Γ
    %  Citation
    %  --------
    %  UQMaT: see CITATION.cff.
    methods (Static)
        function n=idx(ix,iy,Ny), n=(ix-1)*Ny+iy; end
        function [ix,iy]=inv_idx(n,Ny), ix=ceil(n/Ny); iy=n-(ix-1)*Ny; end
        function [kpath,klabels,kticks]=kpath_square(Nk)
            if nargin<1, Nk=100; end
            Ns=ceil(Nk/3);
            G=[0,0];X=[pi,0];M=[pi,pi];
            s1=[linspace(G(1),X(1),Ns)',linspace(G(2),X(2),Ns)'];
            s2=[linspace(X(1),M(1),Ns)',linspace(X(2),M(2),Ns)'];
            s3=[linspace(M(1),G(1),Ns)',linspace(M(2),G(2),Ns)'];
            kpath=[s1;s2(2:end,:);s3(2:end,:)];
            klabels={'\Gamma','X','M','\Gamma'};
            kticks=[1,Ns,2*Ns-1,size(kpath,1)];
        end
        function info(Nx,Ny,bcx,bcy)
            fprintf('\n  Square2D %d×%d=%d sites  bcx=%s bcy=%s\n',...
                Nx,Ny,Nx*Ny,upper(bcx),upper(bcy));
            fprintf('  n(ix,iy) = (ix-1)*Ny + iy  (column-major)\n\n');
        end
    end
end
