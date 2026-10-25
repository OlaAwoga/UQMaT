classdef Junction1D < UQMaT.Models.Base
    % UQMaT.Models.Heterostructure.Junction1D
    %   Two 1D chains (left L and right R) coupled at an interface.
    %
    %   Layout:   [L chain: N_L sites] | Γ | [R chain: N_R sites]
    %   Total:    N = N_L + N_R sites
    %
    %   Block structure:
    %     H = [ H_L,    T_LR ]
    %         [ T_LR†,  H_R  ]
    %
    %   Coupling matrix T_LR:
    %     - Point contact (default): Γ connects site N_L to site N_L+1
    %     - Interfacial spread: Γ connects the last n_contact sites of L
    %       to the first n_contact sites of R (controlled by Gamma_spread)
    %
    %   Each side has independent parameters (system_L, system_R structs).
    %   Anderson disorder on Γ is supported via Gamma_disorder.
    %
    %   Params fields used:
    %     hetero.L : struct  parameters for left chain  (same fields as Params)
    %     hetero.R : struct  parameters for right chain (same fields as Params)
    %     hetero.Gamma        : coupling amplitude (scalar)
    %     hetero.Gamma_spread : number of interface sites coupled (default 1)
    %     hetero.Gamma_disorder: Anderson disorder on Gamma (default 0)
    %     hetero.type         : 'junction1d'
    %
    %   Usage:
    %     p = UQMaT.Params('system','junction1d', ...
    %         'hetero', struct( ...
    %             'type',    'junction1d', ...
    %             'L', struct('N',40,'t',1,'mu',0,'Delta_s',0.3), ...
    %             'R', struct('N',40,'t',1,'mu',-0.5), ...
    %             'Gamma', 0.5, ...
    %             'Gamma_spread', 1, ...
    %             'Gamma_disorder', 0.1))
    %
    %   Physical use cases:
    %     - NS junction (one SC, one normal)
    %     - SS' junction (two different SCs)
    %     - NIN junction (insulating barrier: small Gamma)
    %     - Andreev reflection, subgap states

    %  Citation
    %  --------
    %  NS junction, Andreev reflection:
    %    Andreev A.F., JETP 19, 1228 (1964).
    %    Blonder G.E., Tinkham M., Klapwijk T.M., PRB 25, 4515 (1982).
    %  Subgap (Yu-Shiba-Rusinov) states:
    %    Yu L., Acta Phys. Sin. 21, 75 (1965).
    %    Shiba H., Prog. Theor. Phys. 40, 435 (1968).
    %    Rusinov A.I., JETP Lett. 9, 85 (1969).
    methods
        function obj = Junction1D(p)
            obj.params = p;
            h   = p.hetero;
            N_L = h.L.N;
            N_R = h.R.N;
            obj.Ndim = N_L + N_R;
            obj.H    = obj.buildH();
        end

        function s = description(obj)
            h = obj.params.hetero;
            s = sprintf('Junction1D N_L=%d|Γ=%.3g|N_R=%d', ...
                h.L.N, h.Gamma, h.R.N);
        end

        function H = buildH(obj)
            p = obj.params;
            h = p.hetero;

            % ---- Build left and right Params objects ----
            pL = UQMaT.Models.Heterostructure.Junction1D.make_side_params(h.L, p);
            pR = UQMaT.Models.Heterostructure.Junction1D.make_side_params(h.R, p);

            % ---- Build individual Hamiltonians ----
            sysL = UQMaT.System(pL);
            sysR = UQMaT.System(pR);
            H_L  = sysL.H;
            H_R  = sysR.H;

            dim_L = size(H_L, 1);
            dim_R = size(H_R, 1);
            N_tot = dim_L + dim_R;

            % ---- Coupling matrix T_LR ----
            Gamma   = h.Gamma;
            spread  = 1;
            if isfield(h,'Gamma_spread'), spread = max(1, h.Gamma_spread); end
            Gam_dis = 0;
            if isfield(h,'Gamma_disorder'), Gam_dis = h.Gamma_disorder; end

            T_LR = zeros(dim_L, dim_R);

            % Determine BdG multiplier: if either side is SC, coupling
            % enters in spin-up block (assumes aligned spin basis)
            is_sc = pL.sc || pR.sc;
            is_sp = pL.spin || pR.spin;

            % Site indices at the interface
            % L: last 'spread' electron sites; R: first 'spread' electron sites
            if is_sc
                % BdG: electron block is rows 1:N, spin-up is 1:N/2
                NL_e = pL.N;  % number of electron sites per spin
                NR_e = pR.N;
                block_factor = 4;   % 4N×4N
            elseif is_sp
                NL_e = pL.N;
                NR_e = pR.N;
                block_factor = 2;
            else
                NL_e = pL.N;
                NR_e = pR.N;
                block_factor = 1;
            end

            % Apply Gamma with optional Anderson disorder
            for k = 1:spread
                iL = NL_e - spread + k;   % site in L (counting from interface)
                iR = k;                    % site in R (counting from interface)
                if iL >= 1 && iL <= NL_e && iR >= 1 && iR <= NR_e
                    g = Gamma;
                    if Gam_dis > 0
                        g = g + Gam_dis * (rand - 0.5);
                    end
                    if block_factor == 1
                        T_LR(iL, iR) = -g;
                    elseif block_factor == 2
                        % Couple spin-up and spin-down blocks
                        T_LR(iL,        iR)        = -g;
                        T_LR(iL+NL_e,   iR+NR_e)  = -g;
                    else
                        % BdG: couple electron blocks (spin-up and spin-down)
                        T_LR(iL,            iR)            = -g;
                        T_LR(iL+NL_e,       iR+NR_e)       = -g;
                        % Hole block (PH symmetry): automatic from BdG assembly
                        % but we add it for self-consistency
                        T_LR(iL+2*NL_e,     iR+2*NR_e)     = +conj(g);
                        T_LR(iL+3*NL_e,     iR+3*NR_e)     = +conj(g);
                    end
                end
            end

            % ---- Assemble ----
            H = [ H_L,     T_LR;
                  T_LR',   H_R  ];
        end
    end

    methods (Static)
        function p_side = make_side_params(side_struct, parent_params)
            % Build a Params object for one side of the junction.
            % Inherits unset fields from parent_params.
            f   = fieldnames(side_struct);
            args = {};
            for k = 1:numel(f)
                args{end+1} = f{k}; %#ok<AGROW>
                args{end+1} = side_struct.(f{k}); %#ok<AGROW>
            end
            % Default system to 'chain' if not specified in the struct.
            % For 2D substrates (Chain1D_Square2D), pass system='square'
            % inside the struct to override this default.
            if ~isfield(side_struct,'system')
                args = [{'system','chain'}, args];
            end
            % Inherit eta, EnergyGrid, seed from parent
            for field = {'eta','EnergyGrid','seed'}
                if ~isfield(side_struct, field{1})
                    args{end+1} = field{1}; %#ok<AGROW>
                    args{end+1} = parent_params.(field{1}); %#ok<AGROW>
                end
            end
            p_side = UQMaT.Params(args{:});
        end
    end
end
