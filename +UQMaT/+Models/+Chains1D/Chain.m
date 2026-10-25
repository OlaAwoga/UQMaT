classdef Chain < UQMaT.Models.Base
    % UQMaT.Models.Chains1D.Chain
    %   Uniform 1D tight-binding chain.
    %
    %   H sizes — determined automatically from parameters:
    %     sc=false, spin=false  →  N × N          normal spinless
    %     sc=false, spin=true   →  2N × 2N         spinful normal
    %     sc=true               →  4N × 4N         BdG Ψ=(c↑,c↓,c↑†,c↓†)ᵀ
    %
    %   Superconductivity (sc=true) is activated automatically when any
    %   pairing amplitude is nonzero:
    %     Delta_s ≠ 0  →  singlet s-wave
    %     dvec   ≠ 0  →  triplet p-wave
    %     both         →  mixed
    %
    %   Error is raised if sc=true but all pairing amplitudes are zero.

    %  Citation
    %  --------
    %  Tight-binding chain: Bloch F., Z. Phys. 52, 555 (1929).
    %  Anderson P.W., Phys. Rev. 109, 1492 (1958). [Anderson localisation]
    %  Datta S., "Electronic Transport in Mesoscopic Systems", Cambridge (1995).
    methods
        function obj = Chain(p)
            obj.params = p;
            obj.Ndim   = p.N;
            obj.H      = obj.buildH();
        end

        function s = description(obj)
            p = obj.params;
            if p.sc
                s = sprintf('Chain 1D SC BdG 4×%d | Δ_s=%.3g', ...
                    p.N, p.sc_Delta_s);
            elseif p.spin
                s = sprintf('Chain 1D spinful 2×%d | Bz=%.3g λ_R=%.3g', ...
                    p.N, norm(p.B), p.lambda_R);
            else
                s = sprintf('Chain 1D normal %d | t=%.3g μ=%.3g', ...
                    p.N, p.t, p.mu);
            end
        end

        function H = buildH(obj)
            p   = obj.params;
            N   = p.N;
            bdg = UQMaT.SC.BdG;

            % Validate: sc=true requires at least one nonzero amplitude
            if p.sc && all(p.sc_Delta_s == 0) && p.sc_Delta_t == 0
                error('UQMaT:Chain:noAmplitude', ...
                    ['sc=true but all pairing amplitudes are zero.\n' ...
                     'Set SC.Delta_s or SC.Delta_t to a nonzero value.\n' ...
                     'Example: UQMaT.Params(''system'',''chain'',''SC'',struct(''Delta_s'',0.3))']);
            end

            % Order parameter contributions (all channels)
            [D_s_ops, D_t_ops, h_AM, h_CDW] = UQMaT.OrderParameterSymmetry.OPS.build(p);

            % Scalar kinetic + disorder + CDW (N×N)
            h0 = UQMaT.Models.Chains1D.Chain.hop1D(N, p.t, p.t2, p.bc) ...
                 - p.mu * eye(N) ...
                 + p.onsite_val * eye(N) ...
                 + diag(UQMaT.Disorder.Disorder.onsite(p, N)) ...
                 + UQMaT.Disorder.Disorder.hopping(p, N, p.bc) ...
                 + h_CDW;

            if p.sc
                % ---- BdG: 4N×4N ----
                H_N = kron(eye(2), h0);
                if any(p.B ~= 0)
                    H_N = H_N + bdg.zeeman(p.B, N);
                end
                if p.lambda_R ~= 0
                    A_R = UQMaT.SC.BdG.rashba_hop(N, p.bc);
                    H_N = H_N + bdg.rashba(p.lambda_R, A_R);
                end
                H_N = H_N + h_AM;
                % Combine singlet + triplet pairing blocks
                D = D_s_ops + D_t_ops;
                H = bdg.assemble(H_N, D);

            elseif p.spin
                % ---- Spinful normal: 2N×2N ----
                H = kron(eye(2), h0);
                if any(p.B ~= 0)
                    H = H + bdg.zeeman(p.B, N);
                end
                if p.lambda_R ~= 0
                    A_R = UQMaT.SC.BdG.rashba_hop(N, p.bc);
                    H  = H + bdg.rashba(p.lambda_R, A_R);
                end
                H = H + h_AM;

            else
                % ---- Normal spinless: N×N ----
                H = h0;
            end
        end
    end


    methods (Static)
        function T = hop1D(N, t, t2, bc)
            % Build N×N 1D hopping matrix for given boundary conditions.
            e1 = -t * ones(N-1, 1);
            T  = diag(e1, 1) + diag(e1, -1);
            if t2 ~= 0
                e2 = -t2 * ones(N-2, 1);
                T  = T + diag(e2, 2) + diag(e2, -2);
            end
            if strcmp(bc, 'pbc')
                T(1,N) = -t;   T(N,1) = -t;
                if t2 ~= 0 && N > 3
                    T(1,N-1)=-t2; T(2,N)=-t2;
                    T(N-1,1)=-t2; T(N,2)=-t2;
                end
            end
        end
    end
end
