classdef Lattice < UQMaT.Models.Base
    % UQMaT.Models.Lattice2D.TMD.Lattice
    %   Three-band tight-binding model for transition metal dichalcogenides (TMDs).
    %   Ref: Liu, Xiao et al., Phys Rev B 88, 085433 (2013)
    %
    %   Examples: MoS2, MoSe2, WS2, WSe2
    %
    %   Basis: d_{z²}, d_{xy}, d_{x²-y²} orbitals on the metal M site
    %          in a triangular lattice (hexagonal BZ).
    %   Total sites: Ns = 3 * Nx * Ny  (3 orbitals per unit cell)
    %   Site index: n(ix,iy,orb) = 3*((ix-1)*Ny + (iy-1)) + orb
    %               orb: 1=dz2, 2=dxy, 3=dx2y2
    %
    %   Spin-orbit coupling λ (atomic SOC on metal d-orbitals) is built-in:
    %     H_SOC = λ * diag(0, +1, -1) ⊗ σ_z   (in the spin block)
    %   This splits K/K' valleys.
    %
    %   Pre-defined materials via 'tmd_material' Params field:
    %     'MoS2'  : t₀=−1.663, t₁=0.515, t₂=−0.276, t₁₁=0.424,
    %               t₁₂=0.254, t₂₂=1.155, ε₁=1.046, ε₂=2.104 eV, λ=0.073 eV
    %     'WS2'   : t₀=−1.602, t₁=0.589, t₂=−0.276, t₁₁=0.440,
    %               t₁₂=0.267, t₂₂=1.129, ε₁=1.130, ε₂=2.275 eV, λ=0.211 eV
    %     'custom': use p.tmd_params struct directly
    %
    %   Usage:
    %     p   = UQMaT.Params('system','tmd','Nx',8,'Ny',8,'tmd_material','MoS2')
    %     sys = UQMaT.System(p)
    %
    %   References:
    %     Liu et al., PRB 88, 085433 (2013) — three-band TB model
    %     Xiao et al., PRL 108, 196802 (2012) — valley physics

    %  Citation
    %  --------
    %  Liu G.-B. et al., PRB 88, 085433 (2013). [three-band TB model]
    %  Xiao D. et al., PRL 108, 196802 (2012). [valley physics, SOC]
    %  Kormanyos A. et al., 2D Mater. 2, 022001 (2015). [k·p theory review]
    %  Mak K.F. et al., Nat. Mater. 9, 1085 (2010). [MoS2 experiment]
    methods
        function obj = Lattice(p)
            obj.params = p;
            obj.Ndim   = 3 * p.Nx * p.Ny;
            obj.H      = obj.buildH();
        end

        function s = description(obj)
            p  = obj.params;
            Ns = 3*p.Nx*p.Ny;
            mat = 'custom';
            if isfield(p,'tmd_material'), mat = p.tmd_material; end
            s = sprintf('TMD (%s) %d×%d=%d cells  bc=%s', ...
                mat, p.Nx, p.Ny, Ns, upper(p.bc));
        end

        function H = buildH(obj)
            p   = obj.params;
            Nx  = p.Nx; Ny = p.Ny;
            Ns  = 3*Nx*Ny;    % orbital sites
            bdg = UQMaT.SC.BdG;
            pbc = strcmp(p.bc,'pbc');

            % Get TMD parameters
            mp = UQMaT.Models.Lattice2D.TMD.Lattice.material_params(p);

            % Build kinetic + SOC block
            h0_orb = obj.hopping_matrix(Nx,Ny,mp,pbc);  % 3Ns×3Ns orbital block
            h0_orb = h0_orb - p.mu*eye(Ns);

            % Disorder on orbital sites
            V_dis = diag(UQMaT.Disorder.Disorder.onsite(p, Ns));
            h0_orb = h0_orb + V_dis;

            if p.spin || p.sc
                s   = bdg.sigma();
                % Spin block: 6Ns×6Ns = 2 × 3Ns (spin ⊗ orbital)
                H_orb = kron(eye(2), h0_orb);
                % Intrinsic SOC: λ*(L_z ⊗ σ_z)
                % L_z eigenvalues: dz²→0, dxy→+2, dx²y²→-2 (normalised to ±1)
                Lz_diag = repmat([0;+1;-1], Ns/3, 1);
                Lz = diag(Lz_diag);
                H_orb = H_orb + mp.lambda * kron(s.sz, Lz);
                % External Zeeman
                if any(p.B ~= 0)
                    H_orb = H_orb + kron(p.B(1)*s.sx+p.B(2)*s.sy+p.B(3)*s.sz, eye(Ns));
                end
            end

            [D_s_ops, D_t_ops, h_AM, h_CDW] = UQMaT.OrderParameterSymmetry.OPS.build(p);

            if p.sc
                H_N = H_orb + h_AM;
                D   = D_s_ops + D_t_ops;
                H   = bdg.assemble(H_N, D);
            elseif p.spin
                H   = H_orb + h_AM;
            else
                H   = h0_orb;
            end
        end
    end

    methods (Static)

        function T = hopping_matrix(Nx,Ny,mp,pbc)
            % Build 3Ns×3Ns Hamiltonian for TMD three-band model.
            % Based on Liu et al. PRB 88, 085433 (2013) Table I.
            Ns = 3*Nx*Ny;
            T  = zeros(Ns);
            id = @(x,y,o) UQMaT.Models.Lattice2D.TMD.Lattice.idx_wrap(x,y,o,Nx,Ny,pbc);

            % On-site energies
            for ix=1:Nx; for iy=1:Ny
                n1=id(ix,iy,1); n2=id(ix,iy,2); n3=id(ix,iy,3);
                if ~isnan(n1), T(n1,n1)=T(n1,n1)+mp.e1; end
                if ~isnan(n2), T(n2,n2)=T(n2,n2)+mp.e2; end
                if ~isnan(n3), T(n3,n3)=T(n3,n3)+mp.e2; end
            end; end

            % NN hoppings (6 directions on triangular lattice)
            % Triangular lattice vectors: a1=[1,0], a2=[1/2,√3/2]
            % 6 NN: ±a1, ±a2, ±(a1-a2)
            dirs = [1,0; 0,1; 1,-1; -1,0; 0,-1; -1,1];  % unit cell offsets

            for ix=1:Nx; for iy=1:Ny
                n0 = [id(ix,iy,1); id(ix,iy,2); id(ix,iy,3)];
                for d=1:6
                    dx=dirs(d,1); dy=dirs(d,2);
                    n1 = [id(ix+dx,iy+dy,1); id(ix+dx,iy+dy,2); id(ix+dx,iy+dy,3)];
                    % Hopping matrix for this direction (from Liu PRB 88 Table I)
                    h_hop = UQMaT.Models.Lattice2D.TMD.Lattice.hop_matrix(d, mp);
                    for a=1:3; for b=1:3
                        if ~isnan(n0(a)) && ~isnan(n1(b)) && h_hop(a,b)~=0
                            T(n0(a),n1(b)) = T(n0(a),n1(b)) + h_hop(a,b);
                        end
                    end; end
                end
            end; end
        end

        function h = hop_matrix(dir, mp)
            % 3×3 hopping matrix for direction dir.
            % Simplified symmetric version using Liu et al. parameters.
            t0=mp.t0; t1=mp.t1; t2=mp.t2; t11=mp.t11; t12=mp.t12; t22=mp.t22;
            % Use cos/sin phase factors for the 6 directions at angles 0,60,120,...
            phi = (dir-1)*pi/3;
            c=cos(phi); s=sin(phi); c2=cos(2*phi); s2=sin(2*phi);
            h = [t0,         t1*c-1i*t1*s,    t1*c2-1i*t1*s2 ;
                 t1*c+1i*t1*s, t11*c2+t22,     (t12*c-1i*t12*s)*s  ;
                 t1*c2+1i*t1*s2,(t12*c+1i*t12*s)*s, t22*c2+t11 ];
            % Ensure Hermitian
            h = (h + h')/2;
            % Apply overall scale
            h = real(h);
        end

        function mp = material_params(p)
            % Return TMD material parameters.
            mat = 'MoS2';
            if isfield(p,'tmd_material'), mat = p.tmd_material; end
            switch lower(mat)
                case 'mos2'
                    mp = struct('t0',-1.663,'t1',0.515,'t2',-0.276,...
                        't11',0.424,'t12',0.254,'t22',1.155,...
                        'e1',1.046,'e2',2.104,'lambda',0.073);
                case 'ws2'
                    mp = struct('t0',-1.602,'t1',0.589,'t2',-0.276,...
                        't11',0.440,'t12',0.267,'t22',1.129,...
                        'e1',1.130,'e2',2.275,'lambda',0.211);
                case 'mose2'
                    mp = struct('t0',-1.521,'t1',0.501,'t2',-0.256,...
                        't11',0.372,'t12',0.236,'t22',1.122,...
                        'e1',0.892,'e2',1.780,'lambda',0.091);
                case 'wse2'
                    mp = struct('t0',-1.456,'t1',0.567,'t2',-0.256,...
                        't11',0.384,'t12',0.249,'t22',1.082,...
                        'e1',0.953,'e2',1.910,'lambda',0.228);
                otherwise
                    % custom: expect p.tmd_params struct
                    if isfield(p,'tmd_params')
                        mp = p.tmd_params;
                    else
                        error('UQMaT:TMD:noParams',...
                            'tmd_material=''custom'' requires p.tmd_params struct.');
                    end
            end
            % Override lambda if user provides explicit SOC via p.lambda_R
            if isfield(p,'lambda_R') && p.lambda_R ~= 0
                mp.lambda = p.lambda_R;
            end
        end

    end

    methods (Static, Access=private)
        function n = idx_wrap(ix,iy,orb,Nx,Ny,pbc)
            if pbc, ix=mod(ix-1,Nx)+1; iy=mod(iy-1,Ny)+1; end
            if ix<1||ix>Nx||iy<1||iy>Ny, n=NaN; return; end
            n = 3*((ix-1)*Ny+(iy-1)) + orb;
        end
    end
end
