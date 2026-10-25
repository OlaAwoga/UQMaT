classdef Kitaev < UQMaT.Models.Base
    % UQMaT.Models.Lattice.Kitaev
    %   Kitaev chain — spinless 1D p-wave superconductor.
    %   Always a BdG Hamiltonian; Nambu basis [c_1..c_N, c†_1..c†_N].
    %
    %   Hamiltonian
    %   -----------
    %   H = -t Σ (c†_i c_{i+1} + h.c.) - μ Σ c†_i c_i
    %     + Σ Δ_p (c_i c_{i+1} - c_{i+1} c_i + h.c.)
    %
    %   In Nambu basis (2N × 2N):
    %   H_BdG = [ H_n,  Δ_p ; Δ_p†, -H_n* ]
    %
    %   Topological condition
    %   ---------------------
    %   |μ| < 2t  →  topological  (winding ν=1, Majorana zero modes at OBC edges)
    %   |μ| > 2t  →  trivial      (winding ν=0)
    %   |μ| = 2t  →  phase boundary (gap closes)
    %
    %   With Zeeman field (Bz ≠ 0): 4N × 4N in Nambu ⊗ spin basis.
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params('model','kitaev','N',60,'t',1,'mu',-1,'Delta_p',0.5);
    %     mdl = UQMaT.Models.Lattice.Kitaev(p);
    %     sys = UQMaT.System(mdl);

    methods
        function obj = Kitaev(p)
            if ~isa(p,'UQMaT.Params')
                error('UQMaT:Models:Lattice:Kitaev:badInput', ...
                    'Argument must be a UQMaT.Params object.');
            end
            if ~strcmp(p.model_class,'kitaev')
                error('UQMaT:Models:Lattice:Kitaev:wrongModel', ...
                    'Params.model must be ''kitaev''.');
            end
            obj.params = p;
            obj.Ndim   = obj.ndim();
            obj.H      = obj.buildH();
        end

        function n = ndim(obj)
            n = obj.params.N;
        end

        function s = family(~)
            s = 'lattice';
        end

        function s = description(obj)
            p = obj.params;
            topo = 'trivial';
            if abs(p.mu) < 2*p.t, topo = 'topological'; end
            if abs(p.mu) == 2*p.t, topo = 'phase boundary'; end
            s = sprintf('Kitaev chain, N=%d, bc=%s, t=%.3g μ=%.3g Δ_p=%.3g [%s]', ...
                p.N, upper(p.bc), p.t, p.mu, mean(p.sc_Delta_s_p_vec), topo);
        end

        function H = buildH(obj)
            p  = obj.params;
            N  = p.N;
            Hn = obj.hopMatrix1D(N, p.t, p.bc) - p.mu * eye(N);
            Hn = Hn + obj.disorderMatrix(N);

            % p-wave pairing (antisymmetric)
            Dp = zeros(N);
            for n = 1:N-1
                Dp(n,n+1) =  p.sc_Delta_s_p_vec(n);
                Dp(n+1,n) = -p.sc_Delta_s_p_vec(n);
            end
            if strcmp(p.bc,'pbc')
                Dp(N,1) =  p.sc_Delta_s_p_vec(end);
                Dp(1,N) = -p.sc_Delta_s_p_vec(end);
            end

            if ~p.spin
                H = [Hn, Dp; Dp', -conj(Hn)];
            else
                % Zeeman extension: 4N × 4N
                sz = [1 0;0 -1]; sx = [0 1;1 0]; sy = [0 -1i;1i 0];
                I2 = eye(2);
                Hkin = kron(I2, Hn) + ...
                       kron(p.B(3)*sz + p.B(1)*sx + p.B(2)*sy, eye(N));
                D2N  = kron(1i*sy, Dp);
                H    = [Hkin, D2N; D2N', -conj(Hkin)];
            end
        end
    end

    methods (Access = private)
        function T = hopMatrix1D(~, N, t, bc)
            e1 = -t * ones(N-1,1);
            T  = diag(e1,1) + diag(e1,-1);
            if strcmp(bc,'pbc'), T(1,N) = -t; T(N,1) = -t; end
        end

        function Vd = disorderMatrix(obj, N)
            p = obj.params;
            if p.seed > 0, rng(p.seed); end
            Vd = zeros(N);
            switch p.disorder
                case 'none'
                case 'anderson'
                    Vd = diag(p.W * (rand(N,1) - 0.5));
                case 'binary'
                    Vd = diag(p.W * (2*round(rand(N,1)) - 1) / 2);
                case 'impurity'
                    for k = 1:numel(p.imp_sites)
                        i = p.imp_sites(k);
                        if i >= 1 && i <= N
                            Vd(i,i) = p.imp_potential(k);
                        end
                    end
            end
        end
    end
end
