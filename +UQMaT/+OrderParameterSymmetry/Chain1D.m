classdef Chain1D
    % UQMaT.OrderParameterSymmetry.Chain1D
    %   Order parameter bond structures for the 1D chain.
    %
    %   Each symmetry defines how the order parameter couples to bonds:
    %
    %     onsite      : Δ(i) — on-site, k-independent (s-wave SC)
    %     px          : Δ·[+1,-1]/√2 on bonds (i,i+1) — p_x wave
    %     extended_s  : Δ·[+1,+1]/√2 on bonds (i,i+1) — extended s-wave
    %     dx2         : same as extended_s for 1D (cos k)
    %
    %   Usage:
    %     sym = UQMaT.OrderParameterSymmetry.Chain1D.get('px', N, bc)
    %     sym.bonds    : N×Nbonds  (site, bond) index pairs
    %     sym.weights  : 1×Nbonds complex weights
    %     sym.D        : N×N bond matrix with weights applied
    %
    %   For SC:    D_pairing = Delta * sym.D
    %   For CDW:   h0 += Delta * diag(sym.site_weights)
    %   For AM:    h0 += Delta * diag(sym.sublattice)

    %  Citation
    %  --------
    %  Bond-order SC symmetries:
    %    Sigrist M., Ueda K., Rev. Mod. Phys. 63, 239 (1991). [gap symmetry]
    %  p-wave SC in 1D: Kitaev A.Yu., Phys. Uspekhi 44, 131 (2001).
    methods (Static)

        function sym = get(label, N, bc)
            % Return symmetry struct for the given label.
            if nargin < 3, bc = 'obc'; end
            label = lower(strtrim(label));
            switch label
                case 'onsite'
                    sym = UQMaT.OrderParameterSymmetry.Chain1D.onsite_(N);
                case 'px'
                    sym = UQMaT.OrderParameterSymmetry.Chain1D.px_(N, bc);
                case 'extended_s'
                    sym = UQMaT.OrderParameterSymmetry.Chain1D.extended_s_(N, bc);
                case 'dx2'
                    sym = UQMaT.OrderParameterSymmetry.Chain1D.extended_s_(N, bc);
                otherwise
                    error('UQMaT:OPS:Chain1D:unknown', ...
                        'Unknown symmetry "%s".\nAvailable: onsite, px, extended_s, dx2', label);
            end
            sym.label  = label;
            sym.system = 'chain';
        end

        function list = available()
            list = {'onsite','px','extended_s','dx2'};
        end

    end

    methods (Static, Access = private)

        function sym = onsite_(N)
            % On-site: uniform Δ(i) = 1 on each site
            % SC: s-wave singlet  D_ij = Δ δ_ij
            sym.site_weights = ones(N,1);
            sym.sublattice   = ones(N,1);
            % Bond matrix: diagonal (on-site coupling)
            sym.D = eye(N);
        end

        function sym = px_(N, bc)
            % p_x wave: antisymmetric on bonds (i,i+1)
            % weights: [+1,-1]/√2 for forward/backward bond
            % D(i,i+1) = +1/√2,  D(i+1,i) = -1/√2  (antisymmetric → triplet)
            w = 1/sqrt(2);
            A = zeros(N);
            for i = 1:N-1
                A(i,   i+1) = +w;
                A(i+1, i  ) = -w;
            end
            if strcmp(bc,'pbc')
                A(N, 1) = +w;
                A(1, N) = -w;
            end
            sym.D            = A;
            sym.site_weights = zeros(N,1);   % no on-site component
            sym.sublattice   = (-1).^(1:N)'; % staggered (for AM)
        end

        function sym = extended_s_(N, bc)
            % Extended s-wave: symmetric on bonds (i,i+1)
            % weights: [+1,+1]/√2  (cosine in k)
            % D(i,i+1) = D(i+1,i) = +1/√2  (symmetric → singlet-like)
            w = 1/sqrt(2);
            T = zeros(N);
            for i = 1:N-1
                T(i,   i+1) = w;
                T(i+1, i  ) = w;
            end
            if strcmp(bc,'pbc')
                T(N,1) = w;
                T(1,N) = w;
            end
            sym.D            = T;
            sym.site_weights = zeros(N,1);
            sym.sublattice   = (-1).^(1:N)';
        end

    end
end
