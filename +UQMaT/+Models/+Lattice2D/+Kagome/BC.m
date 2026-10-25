classdef BC
    % UQMaT.Models.Lattice2D.Kagome.BC
    %   Band structure for the kagome lattice via twisted boundary conditions.
    %
    %   BZ high-symmetry path: Γ → M → K → Γ
    %     Γ = (0, 0)
    %     M = (π, 0)
    %     K = (2π/3, 2π/3√3)   (hexagonal BZ corner)
    %
    %   Bloch Hamiltonian H(kx,ky) — 3×3 matrix per k-point (normal state):
    %
    %     H(k) = -t * [ 0,              1+e^{-ikx},     1+e^{-iky}         ]
    %                  [ 1+e^{+ikx},    0,              1+e^{i(kx-ky)}     ]
    %                  [ 1+e^{+iky},    1+e^{-i(kx-ky)}, 0                 ]
    %
    %   Three bands: E_flat = +2t (dispersionless),
    %                E_±(k) = -t ± t√(3 + 2cos(kx) + 2cos(ky) + 2cos(kx-ky))
    %
    %   Usage:
    %     sys = UQMaT.System(UQMaT.Params('system','kagome','Nx',1,'Ny',1,'t',1))
    %     [kd, Ek, kl, kt] = UQMaT.Models.Lattice2D.Kagome.BC.bands(sys)
    %     [kd, Ek, kl, kt] = UQMaT.Models.Lattice2D.Kagome.BC.bands(sys, 300)

    %  Citation
    %  --------
    %  Kagome band structure: Syozi I., Prog. Theor. Phys. 6, 306 (1951).
    %  Flat band: Leykam D. et al., Adv. Phys.: X 3, 1473052 (2018). [review]
    methods (Static)

        function [kdist, Ek, klabels, kticks] = bands(sys, Nk)
            if nargin < 2, Nk = 300; end
            p  = sys.params;
            bc = UQMaT.Models.Lattice2D.Kagome.BC;

            % BZ path: Γ → M → K → Γ
            Ns  = ceil(Nk/3);
            G   = [0,        0       ];
            M   = [pi,       0       ];
            K   = [2*pi/3,   2*pi/3/sqrt(3)];   % hexagonal BZ

            s1 = [linspace(G(1),M(1),Ns)',   linspace(G(2),M(2),Ns)'  ];
            s2 = [linspace(M(1),K(1),Ns)',   linspace(M(2),K(2),Ns)'  ];
            s3 = [linspace(K(1),G(1),Ns)',   linspace(K(2),G(2),Ns)'  ];
            kpts = [s1; s2(2:end,:); s3(2:end,:)];

            klabels = {'\Gamma','M','K','\Gamma'};
            kticks  = [1, Ns, 2*Ns-1, size(kpts,1)];

            Ek = bc.compute_bands(p, kpts);

            % Build kdist
            kdist = zeros(size(kpts,1),1);
            for ik = 2:size(kpts,1)
                kdist(ik) = kdist(ik-1) + norm(kpts(ik,:)-kpts(ik-1,:));
            end
        end

    end

    methods (Static, Access = private)

        function Ek = compute_bands(p, kpts)
            bdg = UQMaT.SC.BdG;
            s   = bdg.sigma();
            t   = p.t;  mu = p.mu;  t2 = p.t2;

            if p.sc,   Nb = 12;  % 4 × 3 sublattice
            elseif p.spin, Nb = 6;
            else,          Nb = 3;
            end

            Nk = size(kpts,1);
            Ek = zeros(Nk, Nb);

            for ik = 1:Nk
                kx = kpts(ik,1);  ky = kpts(ik,2);
                Hk = UQMaT.Models.Lattice2D.Kagome.BC.bloch_h(kx,ky,t,mu,t2);

                if p.sc
                    HN = kron(eye(2), Hk);
                    if any(p.B ~= 0)
                        HN = HN + kron(p.B(1)*s.sx+p.B(2)*s.sy+p.B(3)*s.sz, eye(3));
                    end
                    if p.lambda_R ~= 0
                        Ak  = UQMaT.Models.Lattice2D.Kagome.BC.rashba_bloch(kx,ky);
                        HN  = HN + p.lambda_R * kron(s.sy, Ak);
                    end
                    % Singlet pairing from SC struct
                    if any(p.sc_Delta_s ~= 0)
                        Dk = mean(p.sc_Delta_s) * UQMaT.Models.Lattice2D.Kagome.BC.pairing_bloch(kx,ky,'onsite');
                        D  = kron(1i*s.sy, Dk);
                    else
                        D  = zeros(6);
                    end
                    Hfull = bdg.assemble(HN, D);
                elseif p.spin
                    Hfull = kron(eye(2), Hk);
                    if any(p.B ~= 0)
                        Hfull = Hfull + kron(p.B(1)*s.sx+p.B(2)*s.sy+p.B(3)*s.sz, eye(3));
                    end
                    if p.lambda_R ~= 0
                        Ak    = UQMaT.Models.Lattice2D.Kagome.BC.rashba_bloch(kx,ky);
                        Hfull = Hfull + p.lambda_R * kron(s.sy, Ak);
                    end
                else
                    Hfull = Hk;
                end

                Ek(ik,:) = sort(real(eig((Hfull+Hfull')/2)));
            end
        end

        function Hk = bloch_h(kx, ky, t, mu, t2)
            % 3×3 kagome Bloch Hamiltonian
            % Basis: (A, B, C) per unit cell
            % NN hoppings with phase factors from lattice vectors
            fAB = -t * (1 + exp(-1i*kx));
            fAC = -t * (1 + exp(-1i*ky));
            fBC = -t * (1 + exp(1i*(kx-ky)));

            Hk = [0,       fAB,      fAC;
                  conj(fAB), 0,      fBC;
                  conj(fAC), conj(fBC), 0];

            % On-site: -mu * I
            Hk = Hk - mu*eye(3);

            % NNN: same-sublattice hoppings → -t2*(phase)
            if t2 ~= 0
                phx = exp(1i*kx); phy = exp(1i*ky);
                nnn_diag = -t2 * (phx + conj(phx) + phy + conj(phy));
                Hk = Hk + nnn_diag * eye(3);
            end
        end

        function Ak = rashba_bloch(kx, ky)
            % Antisymmetric Rashba matrix in k-space (3×3)
            fAB =  1i * sin(kx/2);
            fAC =  1i * sin(ky/2);
            fBC =  1i * sin((kx-ky)/2);
            Ak = [0,        fAB,    fAC;
                  -conj(fAB), 0,    fBC;
                  -conj(fAC),-conj(fBC), 0];
        end

        function Dk = pairing_bloch(kx, ky, sym)
            % 3×3 pairing matrix in k-space
            switch sym
                case 'onsite'
                    Dk = eye(3);
                case 'extended_s'
                    w  = 1/sqrt(3);
                    Dk = w * [0,  1+exp(-1i*kx),  1+exp(-1i*ky);
                               1+exp(1i*kx),  0,  1+exp(1i*(kx-ky));
                               1+exp(1i*ky),  1+exp(-1i*(kx-ky)), 0];
                otherwise
                    Dk = eye(3);
            end
        end

    end
end
