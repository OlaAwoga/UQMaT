classdef HNKitaev < UQMaT.Models.Base
    % UQMaT.Models.Chains1D.HNKitaev
    %   Hatano-Nelson-Kitaev (HNK) chain: non-Hermitian 1D topological SC.
    %
    %   Combines:
    %     (1) Hatano-Nelson non-Hermitian hopping (asymmetric t_L ≠ t_R)
    %     (2) Kitaev spinless p-wave pairing Delta_p
    %
    %   Hopping:
    %     t_R = t - t_nh    (rightward)
    %     t_L = t + t_nh    (leftward)
    %     t_nh = 0  →  Hermitian Kitaev chain (t_R = t_L = t)
    %     t_nh ≠ 0  →  non-Hermitian skin effect under OBC
    %     |t_nh| < t  →  both hoppings same sign (weakly NH)
    %     t_nh = t   →  t_R = 0, one-way transport (maximally NH)
    %
    %   Hamiltonian (normal block, N×N):
    %     H_n(i,i+1) = -(t - t_nh)    [rightward]
    %     H_n(i+1,i) = -(t + t_nh)    [leftward]
    %     H_n(i,i)   = -mu
    %
    %   BdG (2N×2N, Nambu basis Ψ = (c_i, c_i†)ᵀ):
    %     H_BdG = [ H_n,   D_p ]
    %             [-D_p*, -H_n*]
    %     D_p(i,i+1) = +Delta_p,  D_p(i+1,i) = -Delta_p
    %
    %   Non-Hermitian winding number (complex energy plane):
    %     ν_NH = (1/2πi) ∮ d/dk [log det H(k)] dk   (base energy E_B)
    %     ν_NH ≠ 0  →  non-Hermitian topological (edge states under OBC)
    %
    %   Spectrum properties:
    %     - PBC: complex eigenvalues forming loops in the complex plane
    %     - OBC: real eigenvalues (skin effect collapses to real axis for
    %            the Kitaev-like phase, complex for other phases)
    %     - Exceptional points: at g = g_c where two eigenvalues coalesce
    %
    %   Phase diagram (g, mu plane at fixed t, Delta_p):
    %     - Non-Hermitian topological SC:  ν_NH = ±1
    %     - Non-Hermitian trivial:         ν_NH = 0
    %     - Exceptional line separating PBC from skin-effect phase
    %
    %   H sizes:
    %     Delta_p = 0  →  N×N    (non-Hermitian normal)
    %     Delta_p ≠ 0  →  2N×2N  (non-Hermitian BdG)
    %
    %   Params fields:
    %     t         : hopping amplitude               (default 1.0)
    %     mu        : chemical potential              (default 0.0)
    %     Delta_p   : p-wave pairing amplitude        (default 0.0)
    %     t_nh      : non-Hermitian asymmetry         (default 0.0)
    %                 t_nh=0 → Hermitian Kitaev chain
    %     bc        : 'obc' | 'pbc'                  (default 'obc')
    %     W         : Anderson disorder on onsite     (default 0.0)
    %
    %   Usage:
    %     % Non-Hermitian Kitaev chain (topological)
    %     p   = UQMaT.Params('system','hnkitaev','N',60,'t',1,'mu',-1,...
    %                        'Delta_p',0.5,'t_nh',0.3)
    %     sys = UQMaT.System(p)
    %
    %     % Hermitian limit (t_nh=0 → standard Kitaev chain)
    %     p   = UQMaT.Params('system','hnkitaev','N',60,'t',1,'mu',-1,...
    %                        'Delta_p',0.5,'t_nh',0)
    %
    %     % Band structure (complex eigenvalues, PBC)
    %     [kn, Ek] = UQMaT.Models.Chains1D.HNKitaev.bands(sys, 400)
    %
    %     % Non-Hermitian winding number
    %     nu_NH = UQMaT.Models.Chains1D.HNKitaev.nh_winding(p)
    %
    %     % Skin effect: compare OBC vs PBC spectra
    %     [E_obc, E_pbc] = UQMaT.Models.Chains1D.HNKitaev.skin_effect(p)
    %
    %  Citation
    %  --------
    %  Non-Hermitian hopping (Hatano-Nelson model):
    %    Hatano N., Nelson D.R., PRL 77, 570 (1996).
    %    Hatano N., Nelson D.R., PRB 56, 8651 (1997).
    %  Non-Hermitian Kitaev chain / winding number:
    %    Lee T.E., PRL 116, 133903 (2016).
    %    Kawabata K. et al., PRB 98, 085116 (2018).
    %    Yao S., Wang Z., PRL 121, 086803 (2018). [non-Hermitian BBC]
    %  Non-Hermitian skin effect:
    %    Zhang K. et al., PRL 125, 126402 (2020).
    %    Okuma N. et al., PRL 124, 086801 (2020).
    %  Exceptional points in SC:
    %    Kawabata K. et al., Nat. Phys. 15, 387 (2019).
    %  Review of non-Hermitian topology:
    %    Bergholtz E.J. et al., Rev. Mod. Phys. 93, 015005 (2021).

    methods
        function obj = HNKitaev(p)
            obj.params = p;
            obj.Ndim   = p.N;
            obj.H      = obj.buildH();
        end

        function s = description(obj)
            p  = obj.params;
            t_nh = p.t_nh;
            dp   = mean(p.Delta_p_vec);
            if dp == 0
                s = sprintf('HN chain N=%d  t=%.3g  mu=%.3g  t_nh=%.3g', ...
                    p.N, p.t, p.mu, t_nh);
            else
                nu = UQMaT.Models.Chains1D.HNKitaev.nh_winding(p);
                s  = sprintf('HNKitaev N=%d  t_nh=%.3g  Δp=%.3g  ν_NH=%d', ...
                    p.N, t_nh, dp, nu);
            end
        end

        function H = buildH(obj)
            p  = obj.params;
            N  = p.N;
            t    = p.t;
            t_nh = p.t_nh;
            mu   = p.mu;

            % ---- Hatano-Nelson kinetic block (N×N, non-Hermitian) ----
            tR = t - t_nh;   % rightward hopping
            tL = t + t_nh;   % leftward  hopping

            Hn = zeros(N);
            for i = 1:N-1
                Hn(i,   i+1) = -tR;   % forward
                Hn(i+1, i  ) = -tL;   % backward
            end
            if strcmp(p.bc, 'pbc')
                Hn(N, 1) = -tR;
                Hn(1, N) = -tL;
            end

            % Chemical potential + onsite disorder
            Hn = Hn - mu*eye(N) + p.onsite_val*eye(N) ...
                 + diag(UQMaT.Disorder.Disorder.onsite(p, N));

            % ---- Normal state only ----
            dp = mean(p.Delta_p_vec);
            if dp == 0
                H = Hn;
                return
            end

            % ---- Non-Hermitian BdG (2N×2N) ----
            % D_p: antisymmetric pairing on bonds
            Dp = zeros(N);
            for n = 1:N-1
                Dp(n,   n+1) = +p.Delta_p_vec(n);
                Dp(n+1, n  ) = -p.Delta_p_vec(n);
            end
            if strcmp(p.bc, 'pbc')
                Dp(N,1) = +p.Delta_p_vec(end);
                Dp(1,N) = -p.Delta_p_vec(end);
            end

            % BdG: [Hn, Dp; -Dp*, -Hn*]
            % Note: NOT assembled via BdG.assemble (which enforces Hermitian form).
            % The HNK BdG breaks Hermiticity via Hn ≠ Hn†.
            H = [ Hn,        Dp;
                 -conj(Dp), -conj(Hn) ];
        end
    end

    methods (Static)

        % ============================================================
        % Band structure — complex eigenvalues vs k (PBC)
        % Returns complex Ek for non-Hermitian spectrum visualisation
        % ============================================================
        function [kn, Ek] = bands(sys, Nk)
            % Bloch Hamiltonian:
            %   H(k) = -(tR e^{ik} + tL e^{-ik}) - mu   [normal block]
            %   tR = t - t_nh,  tL = t + t_nh
            %   BdG adds Delta_p pairing in k-space
            %
            % Returns:
            %   kn : Nk×1  wavenumbers in [-π, π]
            %   Ek : Nk×Nb complex eigenvalues (Nb=1 normal, Nb=2 BdG)
            if nargin < 2, Nk = 400; end
            p  = sys.params;
            t    = p.t;  t_nh = p.t_nh;  mu = p.mu;
            tR = t - t_nh;  tL = t + t_nh;
            Dp = mean(p.Delta_p_vec);
            kn = linspace(-pi, pi, Nk)';
            Nb = 1 + (Dp ~= 0);   % 1 normal, 2 BdG
            Ek = zeros(Nk, Nb);

            for ik = 1:Nk
                k  = kn(ik);
                ek = -(tR*exp(1i*k) + tL*exp(-1i*k)) - mu;
                if Dp == 0
                    Ek(ik,1) = ek;
                else
                    dp_k = 2i*Dp*sin(k);
                    Hk   = [ek, dp_k; -conj(dp_k), -conj(ek)];
                    ev   = eig(Hk);
                    Ek(ik,:) = ev.';
                end
            end
        end

        % ============================================================
        % Non-Hermitian winding number  ν_NH
        % ============================================================
        function nu = nh_winding(p, E_base)
            % Non-Hermitian winding number w.r.t. base energy E_base.
            %   ν_NH = (1/2πi) ∮ ∂_k log det[H(k) - E_B] dk
            %
            % E_base: reference energy (default 0).
            % ν_NH ≠ 0 → topological under OBC (edge states exist).
            % Complex spectrum forms a loop in complex plane encircling E_base.
            if nargin < 2, E_base = 0; end
            Nk = 1000;
            kv = linspace(-pi, pi, Nk+1); kv = kv(1:end-1);
            t  = p.t;  t_nh = p.t_nh;  mu = p.mu;
            tR = t - t_nh;  tL = t + t_nh;
            Dp = mean(p.Delta_p_vec);

            detH = zeros(1, Nk);
            for ik = 1:Nk
                k  = kv(ik);
                ek = -(tR*exp(1i*k) + tL*exp(-1i*k)) - mu;
                if Dp == 0
                    detH(ik) = ek - E_base;
                else
                    dp_k  = 2i*Dp*sin(k);
                    Hk    = [ek - E_base,  dp_k;
                            -conj(dp_k),  -conj(ek) - E_base];
                    detH(ik) = det(Hk);
                end
            end

            % Winding of det[H(k)] around origin
            phase = unwrap(angle(detH));
            nu    = round((phase(end) - phase(1)) / (2*pi));
        end

        % ============================================================
        % Skin effect: compare OBC vs PBC spectra
        % ============================================================
        function [E_obc, E_pbc] = skin_effect(p)
            % Build both OBC and PBC Hamiltonians, return sorted eigenvalues.
            % Under strong non-Hermiticity (g >> 0):
            %   PBC: eigenvalues form a loop in complex plane
            %   OBC: eigenvalues collapse to real axis (skin effect)
            p_obc = UQMaT.Params('system','hnkitaev','N',p.N,'t',p.t,...
                'mu',p.mu,'Delta_p',mean(p.Delta_p_vec),'t_nh',p.t_nh,'bc','obc');
            p_pbc = UQMaT.Params('system','hnkitaev','N',p.N,'t',p.t,...
                'mu',p.mu,'Delta_p',mean(p.Delta_p_vec),'t_nh',p.t_nh,'bc','pbc');
            sys_obc = UQMaT.System(p_obc);
            sys_pbc = UQMaT.System(p_pbc);
            E_obc   = sort(eig(sys_obc.H), 'ComparisonMethod','real');
            E_pbc   = sort(eig(sys_pbc.H), 'ComparisonMethod','real');
        end

        % ============================================================
        % Exceptional points: find g_c where eigenvalues coalesce
        % ============================================================
        function [g_c, E_c] = exceptional_points(p, g_range)
            % Scan g and find exceptional points (EP) where two
            % eigenvalues merge: |dE/dg| → ∞, ImE = 0 splits.
            %   [g_c, E_c] = HNKitaev.exceptional_points(p, [0, 1])
            if nargin < 2, g_range = [0, 0.9*p.t]; end   % t_nh range, avoid |t_nh|>=t
            Ng   = 200;
            g_v  = linspace(g_range(1), g_range(2), Ng);
            Egap = zeros(Ng, 1);

            for ig = 1:Ng
                p_g = UQMaT.Params('system','hnkitaev','N',p.N,'t',p.t,...
                    'mu',p.mu,'Delta_p',mean(p.Delta_p_vec),'t_nh',g_v(ig),'bc','pbc');
                sys = UQMaT.System(p_g);
                E   = sort(eig(sys.H), 'ComparisonMethod','real');
                % Gap = min |E_i - E_j| over all pairs near zero
                E_near = E(abs(real(E)) < 2);
                if numel(E_near) >= 2
                    pairs = abs(diff(E_near));
                    Egap(ig) = min(pairs);
                else
                    Egap(ig) = Inf;
                end
            end

            % EP: local minima of gap close to zero
            [~, idx] = findpeaks(-Egap, 'MinPeakHeight', -0.1*max(Egap));
            g_c = g_v(idx);
            E_c = zeros(size(g_c));
            for k = 1:numel(g_c)
                p_k = UQMaT.Params('system','hnkitaev','N',p.N,'t',p.t,...
                    'mu',p.mu,'Delta_p',mean(p.Delta_p_vec),'t_nh',g_c(k),'bc','pbc');
                sys = UQMaT.System(p_k);
                E_k = eig(sys.H);
                E_c(k) = E_k(abs(real(E_k)) == min(abs(real(E_k))));
            end
        end

        % ============================================================
        % Phase diagram in (g, mu) plane
        % ============================================================
        function [g_v, mu_v, nu_map] = phase_diagram(p, g_range, mu_range, Npts)
            % Compute NH winding number ν_NH(g, mu).
            %   [g_v, mu_v, nu_map] = HNKitaev.phase_diagram(p, [0,2], [-3,3], 40)
            if nargin < 2, g_range  = [0, 0.9]; end   % t_nh range
            if nargin < 3, mu_range = [-3, 3]; end
            if nargin < 4, Npts     = 40; end

            g_v    = linspace(g_range(1),  g_range(2),  Npts);   % t_nh axis
            mu_v   = linspace(mu_range(1), mu_range(2), Npts);
            nu_map = zeros(Npts);

            for ig = 1:Npts
                for im = 1:Npts
                    p_tmp = UQMaT.Params('system','hnkitaev','N',p.N,...
                        't',p.t,'mu',mu_v(im),...
                        'Delta_p',mean(p.Delta_p_vec),...
                        't_nh',g_v(ig),'bc','pbc');
                    nu_map(im, ig) = UQMaT.Models.Chains1D.HNKitaev.nh_winding(p_tmp);
                end
            end
        end

    end
end
