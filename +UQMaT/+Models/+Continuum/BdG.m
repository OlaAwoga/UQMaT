classdef BdG < UQMaT.Models.Base
    % UQMaT.Models.Continuum.BdG
    %   1D continuum BdG model: free-electron dispersion + s-wave or p-wave
    %   pairing, discretised on a uniform real-space grid.
    %
    %   The Hamiltonian is the Bogoliubov-de Gennes equation in real space:
    %
    %     H_BdG ψ = E ψ,   ψ = (u, v)^T
    %
    %     H_n(x) = -ℏ²/(2m*) d²/dx² - μ(x)  [discretised via finite diff.]
    %     Δ(x)   : pairing amplitude (s-wave: scalar; p-wave: gradient term)
    %
    %   Finite-difference discretisation (2nd order):
    %     -d²/dx² f_i ≈ (2f_i - f_{i-1} - f_{i+1}) / dx²
    %     Effective hopping:  t_eff = ℏ²/(2m* dx²)
    %
    %   Parameters (all via UQMaT.Params, plus continuum-specific below)
    %   ------------------------------------------------------------------
    %   Required Params fields:
    %     N       : number of grid points
    %     bc      : 'obc' | 'pbc'
    %     mu      : chemical potential (can be position-dependent: N-vector)
    %
    %   Continuum-specific (set via UQMaT.Params as extra Name-Value):
    %     dx      : grid spacing (default: 1.0, sets energy scale via t_eff)
    %     m_eff   : effective mass in units of m_e (default: 1.0)
    %     Delta   : s-wave pairing (scalar or N-vector)
    %     Delta_p : p-wave pairing strength (scalar; applied as i*Δ_p d/dx)
    %     Bz,Bx,By: Zeeman field (requires 2-component spinor → 4N×4N)
    %     lambda_R: Rashba SOC
    %
    %   Hamiltonian sizes:
    %     no spin, s-wave  :  2N × 2N   Nambu [e, h]
    %     spin + s-wave    :  4N × 4N   Nambu ⊗ spin
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params('model','continuum_bdg','N',200,'mu',1.0,'Delta',0.2);
    %     mdl = UQMaT.Models.Continuum.BdG(p);
    %     sys = UQMaT.System(mdl);

    properties
        dx    (1,1) double = 1.0
        m_eff (1,1) double = 1.0
    end

    methods
        function obj = BdG(p, varargin)
            if ~isa(p,'UQMaT.Params')
                error('UQMaT:Models:Continuum:BdG:badInput', ...
                    'First argument must be a UQMaT.Params object.');
            end
            obj.params = p;
            % Parse extra continuum options
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
            s = sprintf('Continuum BdG 1D, N=%d pts, dx=%.3g, t_eff=%.3g, μ=%.3g, Δ=%.3g', ...
                p.N, obj.dx, t_eff, mean(p.mu), mean(p.sc_Delta_s));
        end

        function H = buildH(obj)
            p     = obj.params;
            N     = p.N;
            t_eff = 1 / (2 * obj.m_eff * obj.dx^2);   % kinetic energy scale

            % Kinetic energy (finite-difference Laplacian)
            e1 = -t_eff * ones(N-1,1);
            Tn = 2*t_eff*eye(N) + diag(e1,1) + diag(e1,-1);
            if strcmp(p.bc,'pbc')
                Tn(1,N) = -t_eff;  Tn(N,1) = -t_eff;
            end

            % Chemical potential (position-dependent allowed)
            if isscalar(p.mu)
                mu_vec = p.mu * ones(N,1);
            else
                mu_vec = p.mu(:);
            end
            Hn = Tn - diag(mu_vec);

            if ~p.spin
                % Spinless s-wave BdG
                if numel(p.sc_Delta_s) == 1
                    D = p.sc_Delta_s * eye(N);
                else
                    D = diag(p.sc_Delta_s);
                end
                H = [Hn, D; D', -conj(Hn)];
            else
                % Spinful: 4N × 4N
                sz = [1 0;0 -1]; sx = [0 1;1 0]; sy = [0 -1i;1i 0];
                I2 = eye(2);
                Hkin = kron(I2, Hn) + ...
                       kron(p.B(3)*sz + p.B(1)*sx + p.B(2)*sy, eye(N));
                if p.lambda_R ~= 0
                    % Rashba: discretised as iλ_R σ_y d/dx
                    dxM = obj.derivMatrix(N, obj.dx, p.bc);
                    Hkin = Hkin + 1i * p.lambda_R * kron(sy, dxM);
                end
                % s-wave singlet pairing
                if numel(p.sc_Delta_s) == 1
                    Dmat = p.sc_Delta_s * eye(N);
                else
                    Dmat = diag(p.sc_Delta_s);
                end
                D2N = kron(1i*sy, Dmat);
                H   = [Hkin, D2N; D2N', -conj(Hkin)];
            end
        end
    end

    methods (Access = private)
        function D = derivMatrix(~, N, dx, bc)
            % First-order finite difference (antisymmetric)
            e1 = ones(N-1,1) / (2*dx);
            D  = diag(e1,1) - diag(e1,-1);
            if strcmp(bc,'pbc')
                D(1,N) = -1/(2*dx);
                D(N,1) =  1/(2*dx);
            end
        end
    end
end
