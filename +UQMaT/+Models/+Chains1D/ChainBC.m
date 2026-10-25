classdef ChainBC
    % UQMaT.Models.Chains1D.ChainBC
    %   Band structure for 1D chain systems via twisted boundary conditions.
    %
    %   H(k) = H_OBC + (-t)*exp(-ik)*|1><N| + (-t)*exp(+ik)*|N><1|
    %
    %   No separate H(k) is written — the wrap-around bond of H_OBC
    %   gets a phase factor exp(±ik). Scanning k ∈ [-π, π] gives bands.
    %
    %   Usage
    %   -----
    %     sys = UQMaT.System(p);   % any 1D system
    %     [kn, Ek] = UQMaT.Models.Chains1D.ChainBC.bands(sys)
    %     [kn, Ek] = UQMaT.Models.Chains1D.ChainBC.bands(sys, Nk)
    %
    %   Works for: chain, ssh, kitaev.
    %   For N=1 unit cell, gives the exact analytical dispersion.

    %  Citation
    %  --------
    %  Twisted boundary conditions: Kohn W., Phys. Rev. 133, A171 (1964).
    %  Hatsugai Y., Phys. Rev. B 48, 11851 (1993). [bulk-edge correspondence]
    methods (Static)

        function [kn, Ek] = bands(sys, Nk)
            if nargin < 2, Nk = 300; end
            p   = sys.params;
            kn  = linspace(-pi, pi, Nk)';
            bc  = UQMaT.Models.Chains1D.ChainBC;
            bdg = UQMaT.SC.BdG;
            s   = bdg.sigma();

            switch p.system
                case 'chain'
                    Ek = bc.bands_chain(p, kn, s, bdg);
                case 'ssh'
                    Ek = bc.bands_ssh(p, kn);
                case 'kitaev'
                    Ek = bc.bands_kitaev(p, kn);
                case 'hnkitaev'
                    % Use HNKitaev.bands directly — returns complex Ek
                    [~, Ek] = UQMaT.Models.Chains1D.HNKitaev.bands(sys, numel(kn));
                otherwise
                    error('UQMaT:ChainBC:unsupported', ...
                        'bands() not supported for system="%s"', p.system);
            end
        end

    end

    methods (Static, Access = private)

        function Ek = bands_chain(p, kn, s, bdg)
            % ε(k) = -2t cos(k) - 2t2 cos(2k) - mu
            Nk = numel(kn);
            if p.sc, Nb=4; elseif p.spin, Nb=2; else, Nb=1; end
            Ek = zeros(Nk, Nb);
            for ik = 1:Nk
                k  = kn(ik);
                ek = -2*p.t*cos(k) - 2*p.t2*cos(2*k) - p.mu + p.onsite_val;
                if ~p.spin && ~p.sc
                    Ek(ik,1) = ek;
                elseif p.spin && ~p.sc
                    HkN = ek*eye(2) + p.B(1)*s.sx + p.B(2)*s.sy + p.B(3)*s.sz;
                    if p.lambda_R~=0
                        HkN = HkN - p.lambda_R*sin(k)*s.sy;
                    end
                    Ek(ik,:) = sort(real(eig(HkN)))';
                else
                    HkN = ek*eye(2) + p.B(1)*s.sx + p.B(2)*s.sy + p.B(3)*s.sz;
                    if p.lambda_R~=0
                        HkN = HkN - p.lambda_R*sin(k)*s.sy;
                    end
                    Dk  = UQMaT.Models.Chains1D.ChainBC.pair_bloch(p, k);
                    Ek(ik,:) = sort(real(eig(bdg.assemble(HkN,Dk))))';
                end
            end
        end

        function Ek = bands_ssh(p, kn)
            % SSH H(k) = [mu, v+w*exp(-ik); v+w*exp(ik), mu]
            Nk = numel(kn);
            Ek = zeros(Nk, 2);
            for ik = 1:Nk
                k   = kn(ik);
                off = p.v + p.w*exp(-1i*k);
                Hk  = [p.mu, off; conj(off), p.mu];
                Ek(ik,:) = sort(real(eig(Hk)))';
            end
        end

        function Ek = bands_kitaev(p, kn)
            % Kitaev H(k): spinless BdG 2×2 when Delta_p≠0, else scalar
            Nk = numel(kn);
            Dp = mean(p.Delta_p_vec);
            if Dp == 0
                Ek = -2*p.t*cos(kn) - p.mu;
            else
                Ek = zeros(Nk, 2);
                for ik = 1:Nk
                    k  = kn(ik);
                    ek = -2*p.t*cos(k) - p.mu + p.onsite_val;
                    Hk = [ek, 2i*Dp*sin(k); -2i*Dp*sin(k), -ek];
                    Ek(ik,:) = sort(real(eig(Hk)))';
                end
            end
        end

        function Dk = pair_bloch(p, k)
            % 2×2 pairing matrix at wavevector k
            % Singlet: Delta_s * iσ_y  (k-independent, on-site)
            % Triplet: (d·σ)(iσ_y) * sin(k)  (p-wave form factor)
            s = UQMaT.SC.BdG.sigma();
            switch p.pairing_mode
                case 'singlet'
                    Dk = mean(p.Delta_s_vec) * (1i*s.sy);
                case 'triplet'
                    dvec = UQMaT.Models.Chains1D.ChainBC.dvec_from_sym(p.sc_sym_t); dx=dvec(1); dy=dvec(2); dz=dvec(3);
                    dsisy = [dz, dx-1i*dy; dx+1i*dy, -dz];
                    Dk = p.Delta_t_val * dsisy * sin(k);
                case 'mixed'
                    dvec = UQMaT.Models.Chains1D.ChainBC.dvec_from_sym(p.sc_sym_t); dx=dvec(1); dy=dvec(2); dz=dvec(3);
                    dsisy = [dz, dx-1i*dy; dx+1i*dy, -dz];
                    Dk = mean(p.Delta_s_vec)*(1i*s.sy) + ...
                         p.Delta_t_val*dsisy*sin(k);
                otherwise
                    Dk = zeros(2);
            end
        end

        function dvec = dvec_from_sym(sym)
            % Return d-vector [dx,dy,dz] for a given triplet symmetry label.
            % Default is [0,0,1] (p_z, standard Kitaev direction).
            switch lower(strtrim(sym))
                case 'px',       dvec = [1,0,0];
                case 'py',       dvec = [0,1,0];
                case 'pz',       dvec = [0,0,1];
                case 'chiral_p', dvec = [1,1i,0]/sqrt(2);
                otherwise,       dvec = [0,0,1];  % safe default
            end
        end

    end
end
