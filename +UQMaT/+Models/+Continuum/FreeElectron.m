classdef FreeElectron < UQMaT.Models.Base
    % UQMaT.Models.Continuum.FreeElectron
    %   1D free-electron (parabolic band) model discretised on a real-space grid.
    %   No pairing; normal state only. Supports Zeeman and Rashba.
    %
    %   Dispersion:  ε(k) = ℏ²k²/(2m*) - μ
    %   Discretised: H_ij = 2t_eff δ_ij - t_eff δ_{i,j±1}
    %                t_eff = ℏ²/(2m* dx²)
    %
    %   Hamiltonian sizes:
    %     no spin  :  N × N
    %     spin     :  2N × 2N
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params('model','continuum_free','N',200,'mu',1.0);
    %     mdl = UQMaT.Models.Continuum.FreeElectron(p);
    %     sys = UQMaT.System(mdl);

    properties
        dx    (1,1) double = 1.0
        m_eff (1,1) double = 1.0
    end

    methods
        function obj = FreeElectron(p, varargin)
            if ~isa(p,'UQMaT.Params')
                error('UQMaT:Models:Continuum:FreeElectron:badInput', ...
                    'First argument must be a UQMaT.Params object.');
            end
            obj.params = p;
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
            obj.Ndim = obj.ndim();
            obj.H    = obj.buildH();
        end

        function n = ndim(obj)
            n = obj.params.N;
        end

        function s = family(~)
            s = 'continuum';
        end

        function s = description(obj)
            p = obj.params;
            t_eff = 1 / (2 * obj.m_eff * obj.dx^2);
            s = sprintf('Free electron 1D, N=%d pts, dx=%.3g, t_eff=%.3g, μ=%.3g', ...
                p.N, obj.dx, t_eff, p.mu);
        end

        function H = buildH(obj)
            p     = obj.params;
            N     = p.N;
            t_eff = 1 / (2 * obj.m_eff * obj.dx^2);

            e1 = -t_eff * ones(N-1,1);
            Hn = 2*t_eff*eye(N) + diag(e1,1) + diag(e1,-1) - p.mu*eye(N);
            if strcmp(p.bc,'pbc')
                Hn(1,N) = -t_eff;  Hn(N,1) = -t_eff;
            end

            if ~p.spin
                H = Hn;
            else
                sz = [1 0;0 -1]; sx = [0 1;1 0]; sy = [0 -1i;1i 0];
                H = kron(eye(2), Hn) + ...
                    kron(p.B(3)*sz + p.B(1)*sx + p.B(2)*sy, eye(N));
                if p.lambda_R ~= 0
                    e1d = ones(N-1,1) / (2*obj.dx);
                    Dx  = diag(e1d,1) - diag(e1d,-1);
                    if strcmp(p.bc,'pbc')
                        Dx(1,N) = -1/(2*obj.dx); Dx(N,1) = 1/(2*obj.dx);
                    end
                    H = H + 1i * p.lambda_R * kron(sy, Dx);
                end
            end
        end
    end
end
