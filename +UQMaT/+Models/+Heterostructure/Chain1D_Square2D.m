classdef Chain1D_Square2D < UQMaT.Models.Base
    % UQMaT.Models.Heterostructure.Chain1D_Square2D
    %   1D chain (or wire array) coupled to a 2D square lattice substrate.
    %
    %   Geometry:
    %     - 2D square substrate: Nx × Ny sites
    %     - 1D chain(s) placed on top, running along x or y
    %     - Each chain site couples to the substrate site directly below it
    %       via coupling Gamma
    %     - Chain position controlled by hetero.chain_row (for x-running chains)
    %       or hetero.chain_col (for y-running chains)
    %
    %   Block structure:
    %     H = [ H_2D,   T_2D_1D ]
    %         [ T_2D_1D†, H_1D  ]
    %
    %   Total sites: Ns_2D + N_chain
    %
    %   hetero struct fields:
    %     chain    : struct  — chain parameters (N,t,mu,Delta_s,...)
    %     substrate: struct  — 2D square lattice parameters (Nx,Ny,t,mu,...)
    %     Gamma         : coupling amplitude (scalar)
    %     Gamma_disorder: Anderson disorder on coupling (default 0)
    %     direction     : 'x' (chain runs along x) | 'y' (along y)
    %     chain_row     : iy position of chain on substrate (for x-direction)
    %     chain_col     : ix position of chain on substrate (for y-direction)
    %     chain_offset  : ix start position of chain (default 1)
    %     chain_length  : number of coupled chain sites (default: full Nx or Ny)
    %
    %   Usage:
    %     p = UQMaT.Params('system','chain1d_square2d', ...
    %         'hetero', struct( ...
    %             'substrate', struct('Nx',20,'Ny',20,'t',1,'mu',0,'Delta_s',0.3), ...
    %             'chain',     struct('N',20,'t',0.8,'mu',-0.5,'B',[0,0,0.4]), ...
    %             'Gamma',     0.3, ...
    %             'direction', 'x', ...
    %             'chain_row', 10, ...
    %             'Gamma_disorder', 0.05))

    %  Citation
    %  --------
    %  Proximity effect: de Gennes P.G., Rev. Mod. Phys. 36, 225 (1964).
    %  Semiconductor-SC proximity: Lutchyn R.M. et al., PRL 105, 077001 (2010).
    %  Oreg Y. et al., PRL 105, 177002 (2010). [wire-on-substrate geometry]
    methods
        function obj = Chain1D_Square2D(p)
            obj.params = p;
            h  = p.hetero;
            % Build substrate Params — inject system='square' into the struct
            % so make_side_params passes it to UQMaT.Params correctly.
            % (pS.system = 'square' after construction silently fails because
            %  system has SetAccess=private on the value class.)
            sub_struct = h.substrate;
            if ~isfield(sub_struct,'system')
                sub_struct.system = 'square';
            end
            pS = UQMaT.Models.Heterostructure.Junction1D.make_side_params(sub_struct, p);
            pC = UQMaT.Models.Heterostructure.Junction1D.make_side_params(h.chain, p);
            Ns_2d  = pS.Nx * pS.Ny;
            obj.Ndim = Ns_2d + pC.N;
            obj.H    = obj.buildH(p, pS, pC);
        end

        function s = description(obj)
            h = obj.params.hetero;
            s = sprintf('Chain1D_Square2D substrate=%dx%d chain_N=%d Γ=%.3g dir=%s', ...
                h.substrate.Nx, h.substrate.Ny, h.chain.N, h.Gamma, h.direction);
        end

        function H = buildH(obj, p, pS, pC)
            h  = p.hetero;
            dir = 'x';
            if isfield(h,'direction'), dir = h.direction; end

            % Build substrate and chain Hamiltonians
            sysS = UQMaT.System(pS);
            sysC = UQMaT.System(pC);
            H_2D = sysS.H;
            H_1D = sysC.H;

            dim_2D = size(H_2D,1);
            dim_1D = size(H_1D,1);

            Ns2d = pS.Nx * pS.Ny;
            Nc   = pC.N;

            % Derive block factors from actual Hamiltonian dimensions.
            % Substrate and chain may be in different bases (e.g. SC substrate
            % → 4×Ns2d, spinful normal chain → 2×Nc) so they must be set
            % independently — never from a shared is_sc / is_sp flag.
            bfS = dim_2D / Ns2d;   % 1 (normal), 2 (spinful), or 4 (BdG)
            bfC = dim_1D / Nc;     % 1 (normal), 2 (spinful), or 4 (BdG)

            % Get coupling parameters
            Gamma  = h.Gamma;
            G_dis  = 0;
            if isfield(h,'Gamma_disorder'), G_dis = h.Gamma_disorder; end

            % Determine chain position on substrate
            chain_row = 1;
            if isfield(h,'chain_row'), chain_row = h.chain_row; end
            chain_col = 1;
            if isfield(h,'chain_col'), chain_col = h.chain_col; end
            chain_offset = 1;
            if isfield(h,'chain_offset'), chain_offset = h.chain_offset; end

            % Build coupling matrix T (dim_2D × dim_1D)
            T = zeros(dim_2D, dim_1D);
            Ny = pS.Ny;

            for k = 1:Nc
                g = Gamma;
                if G_dis > 0, g = g + G_dis*(rand-0.5); end

                if strcmp(dir,'x')
                    % Chain runs along x: chain site k ↔ substrate site (chain_offset+k-1, chain_row)
                    ix = chain_offset + k - 1;
                    iy = chain_row;
                else
                    % Chain runs along y: chain site k ↔ substrate site (chain_col, chain_offset+k-1)
                    ix = chain_col;
                    iy = chain_offset + k - 1;
                end

                if ix >= 1 && ix <= pS.Nx && iy >= 1 && iy <= pS.Ny
                    n2d = (ix-1)*Ny + iy;  % substrate site index
                    n1d = k;               % chain site index
                    % Couple matching spin/particle blocks.
                    % min(bfS,bfC) covers: normal↔normal (1),
                    % spinful↔spinful (2), BdG↔BdG (4), and mixed cases
                    % (e.g. SC substrate bfS=4, spinful chain bfC=2 → 2 blocks).
                    n_blocks = min(bfS, bfC);
                    for b = 0:n_blocks-1
                        row_2d = b * Ns2d + n2d;
                        col_1d = b * Nc   + n1d;
                        T(row_2d, col_1d) = -g;
                    end
                end
            end

            % Assemble
            H = [ H_2D,   T;
                  T',     H_1D ];
        end
    end
end
