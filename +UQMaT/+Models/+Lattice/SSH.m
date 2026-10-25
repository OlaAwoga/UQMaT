classdef SSH < UQMaT.Models.Base
    % UQMaT.Models.Lattice.SSH
    %   Su-Schrieffer-Heeger chain (dimerised 1D lattice).
    %   Two sublattices A and B per unit cell; alternating hoppings v (intra)
    %   and w (inter).
    %
    %   Basis ordering : [A_1, B_1, A_2, B_2, ..., A_N, B_N]  → 2N × 2N
    %
    %   Topological phase
    %   -----------------
    %   w > v  →  topological  (Zak phase = π, edge states at E=0)
    %   w < v  →  trivial      (Zak phase = 0)
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params('model','ssh','N',40,'v',0.3,'w',1.0);
    %     mdl = UQMaT.Models.Lattice.SSH(p);
    %     sys = UQMaT.System(mdl);

    methods
        function obj = SSH(p)
            if ~isa(p,'UQMaT.Params')
                error('UQMaT:Models:Lattice:SSH:badInput', ...
                    'Argument must be a UQMaT.Params object.');
            end
            if ~strcmp(p.model_class,'ssh')
                error('UQMaT:Models:Lattice:SSH:wrongModel', ...
                    'Params.model must be ''ssh''.');
            end
            obj.params = p;
            obj.Ndim   = obj.ndim();
            obj.H      = obj.buildH();
        end

        function n = ndim(obj)
            n = 2 * obj.params.N;   % 2 sites per unit cell
        end

        function s = family(~)
            s = 'lattice';
        end

        function s = description(obj)
            p = obj.params;
            phase = 'trivial';
            if p.w > p.v, phase = 'topological'; end
            s = sprintf('SSH chain, %d cells (%d sites), bc=%s, v=%.3g w=%.3g [%s]', ...
                p.N, 2*p.N, upper(p.bc), p.v, p.w, phase);
        end

        function H = buildH(obj)
            p   = obj.params;
            N   = p.N;
            dim = 2*N;
            H   = zeros(dim);

            for n = 1:N
                iA = 2*n-1; iB = 2*n;
                H(iA,iA) = p.mu;  H(iB,iB) = p.mu;
                H(iA,iB) = -p.v;  H(iB,iA) = -p.v;
                if n < N
                    iA2 = 2*(n+1)-1;
                    H(iB,iA2) = -p.w;  H(iA2,iB) = -p.w;
                end
            end
            if strcmp(p.bc,'pbc')
                H(2*N,1) = -p.w;  H(1,2*N) = -p.w;
            end

            % On-site disorder (applied site-by-site over 2N sites)
            if ~strcmp(p.disorder,'none')
                if p.seed > 0, rng(p.seed); end
                switch p.disorder
                    case 'anderson'
                        v_dis = p.W * (rand(dim,1) - 0.5);
                    case 'binary'
                        v_dis = p.W * (2*round(rand(dim,1)) - 1) / 2;
                    otherwise
                        v_dis = zeros(dim,1);
                end
                H = H + diag(v_dis);
            end
        end
    end
end
