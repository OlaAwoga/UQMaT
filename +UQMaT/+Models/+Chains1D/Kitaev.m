classdef Kitaev < UQMaT.Models.Base
    % UQMaT.Models.Chains1D.Kitaev
    %   Kitaev spinless p-wave superconductor chain.
    %   When Delta_p = 0: normal state N×N.
    %   When Delta_p ≠ 0: spinless BdG 2N×2N, Nambu basis (c, c†).
    %
    %   Topological: |mu| < 2t  (winding ν=1, Majorana zero modes)
    %   Trivial:     |mu| > 2t  (winding ν=0)

    %  Citation
    %  --------
    %  Kitaev A.Yu., Phys. Uspekhi 44, 131 (2001). [Kitaev chain, Majorana modes]
    %  Alicea J., Rep. Prog. Phys. 75, 076501 (2012). [review, topological SC]
    %  Mourik V. et al., Science 336, 1003 (2012). [experimental signature]
    methods
        function obj = Kitaev(p)
            obj.params = p;
            obj.Ndim   = p.N;
            obj.H      = obj.buildH();
        end

        function s = description(obj)
            p  = obj.params;
            dp = mean(p.Delta_p_vec);
            if dp == 0
                s = sprintf('Kitaev normal %d sites', p.N);
            elseif abs(p.mu) < 2*p.t
                s = sprintf('Kitaev BdG 2×%d topological |mu|<2t', p.N);
            elseif abs(p.mu) > 2*p.t
                s = sprintf('Kitaev BdG 2×%d trivial |mu|>2t', p.N);
            else
                s = sprintf('Kitaev BdG 2×%d phase boundary |mu|=2t', p.N);
            end
        end

        function H = buildH(obj)
            p   = obj.params;
            N   = p.N;
            Hn  = UQMaT.Models.Chains1D.Chain.hop1D(N, p.t, 0, p.bc) ...
                  - p.mu*eye(N) ...
                  + p.onsite_val * eye(N) ...
                  + diag(UQMaT.Disorder.Disorder.onsite(p, N)) ...
                  + UQMaT.Disorder.Disorder.hopping(p, N, p.bc);

            dp = mean(p.Delta_p_vec);
            if dp == 0
                H = Hn;
                return
            end

            % p-wave pairing: antisymmetric bond matrix
            Dp = zeros(N);
            for n = 1:N-1
                Dp(n,n+1) = +p.Delta_p_vec(n);
                Dp(n+1,n) = -p.Delta_p_vec(n);
            end
            if strcmp(p.bc,'pbc')
                Dp(N,1) = +p.Delta_p_vec(end);
                Dp(1,N) = -p.Delta_p_vec(end);
            end
            H = [Hn, Dp; Dp', -conj(Hn)];
        end
    end

end
