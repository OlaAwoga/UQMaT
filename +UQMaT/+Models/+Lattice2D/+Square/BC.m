classdef BC
    % UQMaT.Models.Lattice2D.Square.BC
    %   Boundary conditions and band structure via twisted BC for 2D square lattice.
    %
    %   Band structure extraction (no separate H(k)):
    %     H(kx,ky) = H_OBC + phase factors on wrap bonds
    %     T_x(kx)(Nx,1) = -t*exp(-i*kx)   (x wrap bond, forward)
    %     T_y(ky)(Ny,1) = -t*exp(-i*ky)   (y wrap bond, forward)
    %
    %   bcx/bcy combinations:
    %     pbc×pbc → full 2D bands: Γ→X→M→Γ
    %     pbc×obc → cylinder x: bands along kx, Ny sub-bands
    %     obc×pbc → cylinder y: bands along ky, Nx sub-bands
    %     obc×obc → error (no good k)
    %  Citation
    %  --------
    %  Twisted BC for bands: Kohn W., Phys. Rev. 133, A171 (1964).
    %  Square lattice BZ path (Γ→X→M→Γ): Bradley C.J., Cracknell A.P.,
    %    "The Mathematical Theory of Symmetry in Solids", Oxford (1972).
    methods (Static)
        function [kdist,Ek,klabels,kticks]=bands(sys,Nk)
            if nargin<2, Nk=200; end
            p=sys.params; bc=UQMaT.Models.Lattice2D.Square.BC;
            pbc_x=strcmp(p.bcx,'pbc'); pbc_y=strcmp(p.bcy,'pbc');
            if pbc_x&&pbc_y
                Ns=ceil(Nk/3); G=[0,0];X=[pi,0];M=[pi,pi];
                s1=[linspace(G(1),X(1),Ns)',linspace(G(2),X(2),Ns)'];
                s2=[linspace(X(1),M(1),Ns)',linspace(X(2),M(2),Ns)'];
                s3=[linspace(M(1),G(1),Ns)',linspace(M(2),G(2),Ns)'];
                kpts=[s1;s2(2:end,:);s3(2:end,:)];
                klabels={'\Gamma','X','M','\Gamma'};
                kticks=[1,Ns,2*Ns-1,size(kpts,1)];
                Ek=bc.bands_xy(p,kpts);
            elseif pbc_x&&~pbc_y
                kpts=[linspace(-pi,pi,Nk)',zeros(Nk,1)];
                klabels={'-\pi','0','\pi'}; kticks=[1,round(Nk/2),Nk];
                Ek=bc.bands_cyl(p,kpts,'x');
            elseif ~pbc_x&&pbc_y
                kpts=[zeros(Nk,1),linspace(-pi,pi,Nk)'];
                klabels={'-\pi','0','\pi'}; kticks=[1,round(Nk/2),Nk];
                Ek=bc.bands_cyl(p,kpts,'y');
            else
                error('UQMaT:BC:noGoodK',...
                    'bcx=obc and bcy=obc: no good k. Use sys.E for spectrum.');
            end
            kdist=zeros(size(kpts,1),1);
            for ik=2:size(kpts,1), kdist(ik)=kdist(ik-1)+norm(kpts(ik,:)-kpts(ik-1,:)); end
        end
    end
    methods (Static, Access=private)
        function Ek=bands_xy(p,kpts)
            % Both periodic: scalar ε(kx,ky) → small Bloch H
            s=UQMaT.SC.BdG.sigma(); bdg=UQMaT.SC.BdG;
            if p.sc, Nb=4; elseif p.spin, Nb=2; else, Nb=1; end
            Nk=size(kpts,1); Ek=zeros(Nk,Nb);
            for ik=1:Nk
                kx=kpts(ik,1); ky=kpts(ik,2);
                ek=-2*p.tx*cos(kx)-2*p.ty*cos(ky)-4*p.t2*cos(kx)*cos(ky)-p.mu+p.onsite_val;
                if ~p.spin&&~p.sc
                    Ek(ik,1)=ek;
                elseif p.spin&&~p.sc
                    HkN=ek*eye(2)+p.B(1)*s.sx + p.B(2)*s.sy + p.B(3)*s.sz;
                    if p.lambda_R~=0
                        HkN=HkN+p.lambda_R*(sin(ky)*s.sx-sin(kx)*s.sy);
                    end
                    Ek(ik,:)=sort(real(eig(HkN)))';
                else
                    HkN=ek*eye(2)+p.B(1)*s.sx + p.B(2)*s.sy + p.B(3)*s.sz;
                    if p.lambda_R~=0
                        HkN=HkN+p.lambda_R*(sin(ky)*s.sx-sin(kx)*s.sy);
                    end
                    Dk=UQMaT.Models.Lattice2D.Square.BC.pair_bloch(p,kx,ky);
                    Hbdg=bdg.assemble(HkN,Dk);
                    Ek(ik,:)=sort(real(eig(Hbdg)))';
                end
            end
        end
        function Ek=bands_cyl(p,kpts,dir)
            % One direction periodic: build real-space H in open direction
            s=UQMaT.SC.BdG.sigma(); bdg=UQMaT.SC.BdG;
            hop=@UQMaT.Models.Chains1D.Chain.hop1D;
            if strcmp(dir,'x'), Nopen=p.Ny; T_open=hop(p.Ny,p.ty,0,'obc');
            else,               Nopen=p.Nx; T_open=hop(p.Nx,p.tx,0,'obc'); end
            if p.sc, Nblk=4*Nopen; elseif p.spin, Nblk=2*Nopen; else, Nblk=Nopen; end
            Nk=size(kpts,1); Ek=zeros(Nk,Nblk);
            for ik=1:Nk
                kx=kpts(ik,1); ky=kpts(ik,2);
                if strcmp(dir,'x'), ek_k=-2*p.tx*cos(kx)-p.mu+p.onsite_val;
                else,               ek_k=-2*p.ty*cos(ky)-p.mu+p.onsite_val; end
                if ~p.spin&&~p.sc
                    H_ky=ek_k*eye(Nopen)+T_open;
                elseif p.spin&&~p.sc
                    H_N=kron(eye(2),ek_k*eye(Nopen)+T_open);
                    if any(p.B ~= 0)
                        H_N=H_N+kron(p.B(1)*s.sx + p.B(2)*s.sy + p.B(3)*s.sz,eye(Nopen));
                    end
                    if p.lambda_R~=0
                        if strcmp(dir,'x')
                            H_N=H_N-p.lambda_R*sin(kx)*kron(s.sy,eye(Nopen))...
                                +p.lambda_R*kron(s.sx,hop(Nopen,p.lambda_R,0,'obc'));
                        else
                            H_N=H_N+p.lambda_R*sin(ky)*kron(s.sx,eye(Nopen))...
                                -p.lambda_R*kron(s.sy,hop(Nopen,p.lambda_R,0,'obc'));
                        end
                    end
                    H_ky=H_N;
                else
                    H_N=kron(eye(2),ek_k*eye(Nopen)+T_open);
                    if any(p.B ~= 0)
                        H_N=H_N+kron(p.B(1)*s.sx + p.B(2)*s.sy + p.B(3)*s.sz,eye(Nopen));
                    end
                    if p.lambda_R~=0
                        if strcmp(dir,'x')
                            H_N=H_N-p.lambda_R*sin(kx)*kron(s.sy,eye(Nopen));
                        else
                            H_N=H_N+p.lambda_R*sin(ky)*kron(s.sx,eye(Nopen));
                        end
                    end
                    Ds=mean(p.Delta_s_vec);
                    Dk_s=Ds*kron(1i*s.sy,eye(Nopen));
                    H_ky=bdg.assemble(H_N,Dk_s);
                end
                Ek(ik,:)=sort(real(eig(H_ky)))';
            end
        end
        function Dk=pair_bloch(p,kx,ky)
            s=UQMaT.SC.BdG.sigma();
            switch p.pairing_mode
                case 'singlet'
                    Dk=mean(p.Delta_s_vec)*(1i*s.sy);
                case 'triplet'
                    dx=0;dy=0;dz=0; % dvec removed — use SC struct
                    dsisy=[dz,dx-1i*dy;dx+1i*dy,-dz];
                    fk=sin(kx)+1i*sin(ky);
                    Dk=p.Delta_t_val*dsisy*fk;
                case 'mixed'
                    dx=0;dy=0;dz=0; % dvec removed — use SC struct
                    dsisy=[dz,dx-1i*dy;dx+1i*dy,-dz];
                    fk=sin(kx)+1i*sin(ky);
                    Dk=mean(p.Delta_s_vec)*(1i*s.sy)+p.Delta_t_val*dsisy*fk;
                otherwise, Dk=zeros(2);
            end
        end
    end
end
