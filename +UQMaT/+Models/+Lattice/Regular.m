classdef Regular < UQMaT.Models.Base
    % UQMaT.Models.Lattice.Regular
    %   Uniform 1D tight-binding chain.
    %
    %   Hamiltonian sizes
    %   -----------------
    %   sc=false, spin=false  :  N × N     normal state
    %   sc=false, spin=true   :  2N × 2N   spinful normal + Zeeman + Rashba
    %   sc=true               :  4N × 4N   BdG in Nambu basis
    %                                       Ψ = (c↑, c↓, c↑†, c↓†)ᵀ
    %
    %   Superconductivity uses the full Nambu basis regardless of pairing type.
    %   Zeeman and Rashba enter the normal block H_N; pairing enters Δ̂.
    %
    %   Usage
    %   -----
    %     % Normal state
    %     p = UQMaT.Params('N',80,'t',1,'mu',0);
    %
    %     % Singlet s-wave SC (Nambu, 4N×4N)
    %     p = UQMaT.Params('N',80,'sc',true,'pairing','singlet','Delta_s',0.3);
    %
    %     % SC + Zeeman + Rashba (topological SC candidate)
    %     p = UQMaT.Params('N',80,'Delta_s',0.3,'Bz',0.5,'lambda_R',0.4);
    %
    %     % Triplet p-wave SC
    %     p = UQMaT.Params('N',80,'pairing','triplet','dvec',[0 0 1],'Delta_t',0.4);
    %
    %     % Mixed singlet + triplet
    %     p = UQMaT.Params('N',80,'pairing','mixed','Delta_s',0.2, ...
    %                      'dvec',[0 0 1],'Delta_t',0.3);
    %
    %     mdl = UQMaT.Models.Lattice.Regular(p);
    %     sys = UQMaT.System(p);    % System routes automatically

    methods
        function obj = Regular(p)
            if ~isa(p,'UQMaT.Params')
                error('UQMaT:Lattice:Regular:badInput', ...
                    'Argument must be a UQMaT.Params object.');
            end
            if ~strcmp(p.model_class,'regular') && ~strcmp(p.system,'chain')
                error('UQMaT:Lattice:Regular:wrongModel', ...
                    'Params.model must be ''regular''.');
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
            if p.bdg
                base = sprintf('Regular SC 4×%d [%s]', p.N, p.pairing_mode_mode);
            elseif p.spin
                base = sprintf('Regular spinful 2×%d', p.N);
            else
                base = sprintf('Regular %d sites', p.N);
            end
            extras = {};
            if p.B(3)~=0||p.B(1)~=0||p.B(2)~=0
                extras{end+1} = sprintf('B=[%.2g,%.2g,%.2g]',p.B(1),p.B(2),p.B(3));
            end
            if p.lambda_R~=0
                extras{end+1} = sprintf('λ_R=%.3g',p.lambda_R);
            end
            if ~strcmp(p.disorder,'none')
                extras{end+1} = sprintf('%s W=%.3g',p.disorder,p.W);
            end
            if isempty(extras)
                s = base;
            else
                s = [base ' | ' strjoin(extras,' ')];
            end
        end

        function H = buildH(obj)
            p  = obj.params;
            N  = p.N;
            nb = UQMaT.Models.Lattice.Nambu.Basis;

            % --- Scalar kinetic + disorder (N×N) ---
            h_k = UQMaT.Models.Lattice.BoundaryConditions.BC.hop1D(N, p.t, p.t2, p.bc) - p.mu * eye(N);
            h_k = h_k + obj.disorderDiag(N);

            if ~p.spin && ~p.bdg
                % ---- Normal spinless (N×N) ----
                H = h_k;

            elseif p.spin && ~p.bdg
                % ---- Spinful normal (2N×2N) ----
                H_N = kron(eye(2), h_k);
                if p.B(1)~=0||p.B(2)~=0||p.B(3)~=0
                    H_N = H_N + nb.zeeman(p.B(1), p.B(2), p.B(3), N);
                end
                if p.lambda_R ~= 0
                    T_R = UQMaT.Models.Lattice.BoundaryConditions.BC.hop1D(N, p.lambda_R, 0, p.bc);
                    H_N = H_N + nb.rashba(p.lambda_R, T_R);
                end
                H = H_N;

            else
                % ---- Superconducting: Nambu basis 4N×4N ----
                % Normal block H_N (2N×2N)
                H_N = kron(eye(2), h_k);
                if p.B(1)~=0||p.B(2)~=0||p.B(3)~=0
                    H_N = H_N + nb.zeeman(p.B(1), p.B(2), p.B(3), N);
                end
                if p.lambda_R ~= 0
                    T_R = UQMaT.Models.Lattice.BoundaryConditions.BC.hop1D(N, p.lambda_R, 0, p.bc);
                    H_N = H_N + nb.rashba(p.lambda_R, T_R);
                end

                % Pairing block Δ̂ (2N×2N)
                Delta_hat = obj.buildPairing(N, p);

                % Assemble BdG
                H = nb.assemble(H_N, Delta_hat);
            end
        end

        function disp(obj)
            p = obj.params;
            fprintf('  UQMaT.Models.Lattice.Regular\n');
            fprintf('  %s\n', obj.description());
            fprintf('  H: %d×%d\n', size(obj.H,1), size(obj.H,2));
            if p.bdg
                nb = UQMaT.Models.Lattice.Nambu.Basis;
                nb.ph_sym(obj.H);
                nb.tr_sym(obj.H);
            end
        end
    end

    % ====================================================================
    methods (Access = private)

        function D = buildPairing(obj, N, p)
            nb = UQMaT.Models.Lattice.Nambu.Basis;
            switch p.pairing_mode_mode
                case 'singlet'
                    D = nb.singlet_pairing(p.sc_Delta_s_s_vec, N);

                case 'triplet'
                    D = nb.triplet_pairing([0,0,0], N, p.sc_Delta_s_t_val, p.bc);

                case 'mixed'
                    D_s = nb.singlet_pairing(p.sc_Delta_s_s_vec, N);
                    D_t = nb.triplet_pairing([0,0,0], N, p.sc_Delta_s_t_val, p.bc);
                    D   = D_s + D_t;
            end
        end

        function Vd = disorderDiag(obj, N)
            p = obj.params;
            if p.seed > 0, rng(p.seed); end
            switch p.disorder
                case 'none'
                    Vd = zeros(N);
                case 'anderson'
                    Vd = diag(p.W * (rand(N,1) - 0.5));
                case 'binary'
                    Vd = diag(p.W * (2*round(rand(N,1)) - 1) / 2);
                case 'impurity'
                    Vd = zeros(N);
                    for k = 1:numel(p.imp_sites)
                        i = p.imp_sites(k);
                        if i >= 1 && i <= N
                            Vd(i,i) = p.imp_potential(k);
                        end
                    end
                otherwise
                    Vd = zeros(N);
            end
        end

    end
end
