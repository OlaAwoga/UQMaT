classdef BC
    % UQMaT.Models.Lattice3D.Cubic.BC
    %   Band structure for the simple cubic lattice.
    %   BZ path: Γ → X → M → R → Γ → M
    %     Γ = (0,0,0)
    %     X = (π,0,0)
    %     M = (π,π,0)
    %     R = (π,π,π)

    %  Citation
    %  --------
    %  Cubic BZ path (Γ→X→M→R): Bradley C.J., Cracknell A.P.,
    %    "Mathematical Theory of Symmetry in Solids", Oxford (1972).
    methods (Static)

        function [kdist, Ek, klabels, kticks] = bands(sys, Nk)
            if nargin < 2, Nk = 300; end
            p  = sys.params;
            tx = p.tx; ty = p.ty; tz = p.tz;
            mu = p.mu;  t2 = p.t2;

            Ns = ceil(Nk/5);
            G  = [0,0,0]; X = [pi,0,0]; M = [pi,pi,0]; R = [pi,pi,pi];
            segs = {G,X; X,M; M,R; R,G; G,M};
            lbl  = {'\Gamma','X','M','R','\Gamma','M'};

            kpts = []; tks = [1];
            for k = 1:5
                a = segs{k,1}; b = segs{k,2};
                seg = [linspace(a(1),b(1),Ns)', linspace(a(2),b(2),Ns)', linspace(a(3),b(3),Ns)'];
                if k>1, seg=seg(2:end,:); end
                kpts = [kpts; seg]; %#ok<AGROW>
                tks(end+1) = size(kpts,1); %#ok<AGROW>
            end
            klabels = lbl;  kticks = tks;

            % Normal dispersion
            Nkk = size(kpts,1);
            if p.sc,   Nb=8; elseif p.spin, Nb=2; else, Nb=1; end
            Ek = zeros(Nkk, Nb);

            bdg = UQMaT.SC.BdG;  s = bdg.sigma();
            for ik = 1:Nkk
                kx=kpts(ik,1); ky=kpts(ik,2); kz=kpts(ik,3);
                ek = -2*tx*cos(kx) - 2*ty*cos(ky) - 2*tz*cos(kz) - mu;
                if p.sc
                    HN = ek*eye(2);
                    if any(p.B~=0), HN=HN+p.B(1)*s.sx+p.B(2)*s.sy+p.B(3)*s.sz; end
                    D  = mean(p.sc_Delta_s)*eye(2);
                    Hf = bdg.assemble(HN, kron(1i*s.sy, D));
                elseif p.spin
                    Hf = ek*eye(2);
                    if any(p.B~=0), Hf=Hf+p.B(1)*s.sx+p.B(2)*s.sy+p.B(3)*s.sz; end
                else
                    Hf = ek;
                end
                Ek(ik,:) = sort(real(eig(Hf)));
            end

            kdist = zeros(Nkk,1);
            for ik=2:Nkk, kdist(ik)=kdist(ik-1)+norm(kpts(ik,:)-kpts(ik-1,:)); end
        end

    end
end
