classdef Disorder
    % UQMaT.Disorder.Disorder
    %   Disorder potential generators for 1D and 2D lattice models.
    %
    %   Public static methods (called from model buildH):
    %     onsite(p, N)                    → N×1 diagonal potential
    %     impurity_2d_or_bulk(p,Ns,Nx,Ny) → Ns×1 for 2D lattice
    %     hopping(p, N, bc)               → N×N off-diagonal perturbation
    %     magnetic(p, N)                  → 2N×2N magnetic spin block (1D)
    %     magnetic_2d(p, Nx, Ny)          → 2Ns×2Ns magnetic spin block (2D)
    %     vacancy(p, N, is2d, Nx, Ny)     → N×1 large-penalty vector
    %     impurity_2d(p, Nx, Ny)          → Ns×1 for 2D impurity geometries
    %
    %   Disorder types (p.disorder):
    %     'anderson'    : ε_i ~ Uniform[-W/2,+W/2], zero mean
    %     'gaussian'    : ε_i ~ N(0,W²), zero mean
    %     'binary'      : ε_i ∈ {-W/2,+W/2}
    %     'hop_disorder': δt_ij random (off-diagonal)
    %     'impurity'    : controlled impurities via imp_pos/potential/magnetic
    %     'vacancy'     : sites at imp_pos removed (large on-site penalty)
    %
    %   Impurity geometry (2D, p.imp_geometry):
    %     'point'   : site(s) at imp_pos
    %     'hline'   : horizontal line at iy = imp_pos(2)
    %     'vline'   : vertical line at ix = imp_pos(1)
    %     'line'    : straight line from imp_pos(1,:) to imp_pos(2,:)
    %     'cluster' : all sites within radius imp_sigma of imp_pos

    %  Citation
    %  --------
    %  Anderson disorder: Anderson P.W., Phys. Rev. 109, 1492 (1958).
    %  Binary alloy: Ziman J.M., "Models of Disorder", Cambridge (1979).
    %  Magnetic impurities (Yu-Shiba-Rusinov):
    %    Yu L., Acta Phys. Sin. 21, 75 (1965).
    %    Shiba H., Prog. Theor. Phys. 40, 435 (1968).
    %    Rusinov A.I., JETP Lett. 9, 85 (1969).
    %  Gaussian spread model: Datta S., "Electronic Transport in Mesoscopic
    %    Systems", Cambridge (1995).
    methods (Static)

        % ----------------------------------------------------------------
        % onsite — main dispatcher for 1D models
        % ----------------------------------------------------------------
        function V = onsite(p, N)
            if p.seed > 0, rng(p.seed); end
            d = UQMaT.Disorder.Disorder;
            switch p.disorder
                case 'none'
                    V = zeros(N, 1);
                case 'anderson'
                    V = d.anderson(N, p.W);
                case 'gaussian'
                    V = d.gaussian(N, p.W);
                case 'binary'
                    V = d.binary(N, p.W);
                case {'impurity','vacancy'}
                    V = d.impurity_onsite(p, N);
                otherwise
                    V = zeros(N, 1);
            end
        end

        % ----------------------------------------------------------------
        % impurity_2d_or_bulk — main dispatcher for 2D lattice
        % ----------------------------------------------------------------
        function V = impurity_2d_or_bulk(p, Ns, Nx, Ny)
            if p.seed > 0, rng(p.seed); end
            d = UQMaT.Disorder.Disorder;
            switch p.disorder
                case 'none'
                    V = zeros(Ns, 1);
                case 'anderson'
                    V = d.anderson(Ns, p.W);
                case 'gaussian'
                    V = d.gaussian(Ns, p.W);
                case 'binary'
                    V = d.binary(Ns, p.W);
                case 'impurity'
                    V = d.impurity_2d(p, Nx, Ny);
                case 'vacancy'
                    V = d.vacancy(p, Ns, true, Nx, Ny);
                otherwise
                    V = zeros(Ns, 1);
            end
        end

        % ----------------------------------------------------------------
        % hopping — off-diagonal bond disorder (N×N matrix)
        % ----------------------------------------------------------------
        function dT = hopping(p, N, bc)
            if ~strcmp(p.disorder, 'hop_disorder')
                dT = zeros(N);
                return
            end
            if p.seed > 0, rng(p.seed); end
            e  = p.W * (rand(N-1,1) - 0.5);
            e  = e - mean(e);
            dT = diag(e,1) + diag(e,-1);
            if strcmp(bc, 'pbc')
                dT(1,N) = p.W*(rand-0.5);
                dT(N,1) = dT(1,N);
            end
        end

        % ----------------------------------------------------------------
        % magnetic — 2N×2N magnetic impurity spin block (1D)
        % ----------------------------------------------------------------
        function HM = magnetic(p, N)
            s  = UQMaT.SC.BdG.sigma();
            VM = UQMaT.Disorder.Disorder.mag_onsite(p, N);
            HM = kron(s.sz, diag(VM));
        end

        % ----------------------------------------------------------------
        % impurity_2d — Ns×1 potential for 2D geometries
        % ----------------------------------------------------------------
        function V = impurity_2d(p, Nx, Ny)
            Ns   = Nx * Ny;
            V    = zeros(Ns, 1);
            if isempty(p.imp_pos), return; end
            d    = UQMaT.Disorder.Disorder;
            geom = lower(strtrim(p.imp_geometry));

            switch geom
                case 'point'
                    pos = d.parse_pos_2d(p.imp_pos);
                    for k = 1:size(pos,1)
                        V = V + d.spread_2d(pos(k,:), p.imp_potential_val(k), ...
                            p.imp_sigma, Nx, Ny, 'gaussian');
                    end

                case 'hline'
                    iy = max(1, min(Ny, round(p.imp_pos(end))));
                    for ix = 1:Nx
                        n = (ix-1)*Ny + iy;
                        V(n) = V(n) + p.imp_potential_val(1);
                    end

                case 'vline'
                    ix = max(1, min(Nx, round(p.imp_pos(1))));
                    for iy = 1:Ny
                        n = (ix-1)*Ny + iy;
                        V(n) = V(n) + p.imp_potential_val(1);
                    end

                case 'line'
                    pos = d.parse_pos_2d(p.imp_pos);
                    if size(pos,1) < 2
                        error('UQMaT:Disorder:line', ...
                            'imp_geometry=''line'' needs imp_pos with 2 rows.');
                    end
                    p1    = pos(1,:);  p2 = pos(2,:);
                    Nstep = max(abs(p2-p1)) + 1;
                    for st = 0:Nstep-1
                        frac = st / max(Nstep-1, 1);
                        ix   = max(1,min(Nx, round(p1(1)+frac*(p2(1)-p1(1)))));
                        iy   = max(1,min(Ny, round(p1(2)+frac*(p2(2)-p1(2)))));
                        n    = (ix-1)*Ny + iy;
                        V(n) = V(n) + p.imp_potential_val(1);
                    end

                case 'cluster'
                    pos = d.parse_pos_2d(p.imp_pos);
                    for k = 1:size(pos,1)
                        V = V + d.spread_2d(pos(k,:), p.imp_potential_val(k), ...
                            p.imp_sigma, Nx, Ny, 'hard');
                    end

                otherwise
                    error('UQMaT:Disorder:badGeometry', ...
                        'imp_geometry must be point|hline|vline|line|cluster. Got "%s"', geom);
            end
        end

        % ----------------------------------------------------------------
        % magnetic_2d — 2Ns×2Ns magnetic impurity spin block (2D)
        % ----------------------------------------------------------------
        function HM = magnetic_2d(p, Nx, Ny)
            s  = UQMaT.SC.BdG.sigma();
            VM = UQMaT.Disorder.Disorder.mag_2d(p, Nx, Ny);
            HM = kron(s.sz, diag(VM));
        end

        % ----------------------------------------------------------------
        % kagome_onsite — disorder for kagome (3*Nx*Ny sites)
        % Thin wrapper: same logic as onsite but Ns = 3*Nx*Ny
        % ----------------------------------------------------------------
        function V = kagome_onsite(p, Ns)
            % Same dispatch as onsite — kagome has no special treatment needed
            V = UQMaT.Disorder.Disorder.onsite(p, Ns);
        end

        % ----------------------------------------------------------------
        % vacancy — large on-site penalty at imp_pos sites
        % ----------------------------------------------------------------
        function V = vacancy(p, N, is2d, Nx, Ny)
            LARGE = 1e6;
            V     = zeros(N, 1);
            if isempty(p.imp_pos), return; end
            d = UQMaT.Disorder.Disorder;
            if is2d
                pos = d.parse_pos_2d(p.imp_pos);
                for k = 1:size(pos,1)
                    ix = max(1,min(Nx,round(pos(k,1))));
                    iy = max(1,min(Ny,round(pos(k,2))));
                    n  = (ix-1)*Ny + iy;
                    V(n) = LARGE;
                end
            else
                for k = 1:numel(p.imp_pos)
                    i = round(p.imp_pos(k));
                    if i >= 1 && i <= N
                        V(i) = LARGE;
                    end
                end
            end
        end

    end   % methods (Static)

    methods (Static, Access = private)

        % ---- Bulk disorder generators --------------------------------
        function V = anderson(N, W)
            v = W*(rand(N,1)-0.5);
            V = v - mean(v);
        end

        function V = gaussian(N, W)
            v = W*randn(N,1);
            V = v - mean(v);
        end

        function V = binary(N, W)
            V = W*(2*round(rand(N,1))-1)/2;
        end

        % ---- 1D impurity with optional Gaussian spread ---------------
        function V = impurity_onsite(p, N)
            V = zeros(N,1);
            if isempty(p.imp_pos), return; end
            d = UQMaT.Disorder.Disorder;
            if strcmp(p.disorder,'vacancy')
                for k = 1:numel(p.imp_pos)
                    i = max(1,min(N,round(p.imp_pos(k))));
                    V(i) = 1e6;
                end
                return
            end
            for k = 1:numel(p.imp_pos)
                i0  = p.imp_pos(k);
                V0  = p.imp_potential_val(k);
                sig = p.imp_sigma;
                if sig <= 0
                    i = max(1,min(N,round(i0)));
                    V(i) = V(i) + V0;
                else
                    x = (1:N)';
                    g = V0 * exp(-(x-i0).^2/(2*sig^2));
                    if sum(g) > 0, g = g*(V0/sum(g)); end
                    V = V + g;
                end
            end
        end

        % ---- Magnetic impurity (1D): N×1 exchange per site ----------
        function VM = mag_onsite(p, N)
            VM = zeros(N,1);
            if isempty(p.imp_pos), return; end
            for k = 1:numel(p.imp_pos)
                i0  = p.imp_pos(k);
                J   = p.imp_magnetic_val(k);
                sig = p.imp_sigma;
                if sig <= 0
                    i = max(1,min(N,round(i0)));
                    VM(i) = VM(i) + J;
                else
                    x  = (1:N)';
                    g  = J * exp(-(x-i0).^2/(2*sig^2));
                    if sum(g) > 0, g = g*(J/sum(g)); end
                    VM = VM + g;
                end
            end
        end

        % ---- 2D Gaussian / hard-disk spread -------------------------
        function V = spread_2d(pos, V0, sigma, Nx, Ny, mode)
            Ns  = Nx*Ny;
            V   = zeros(Ns,1);
            ix0 = pos(1);  iy0 = pos(2);
            for ix = 1:Nx
                for iy = 1:Ny
                    r2 = (ix-ix0)^2 + (iy-iy0)^2;
                    n  = (ix-1)*Ny + iy;
                    if strcmp(mode,'hard')
                        if r2 <= sigma^2
                            V(n) = V(n) + V0;
                        end
                    else
                        if sigma <= 0
                            if ix==round(ix0) && iy==round(iy0)
                                V(n) = V(n) + V0;
                            end
                        else
                            V(n) = V(n) + V0*exp(-r2/(2*sigma^2));
                        end
                    end
                end
            end
            if sigma > 0 && strcmp(mode,'gaussian')
                s = sum(V);
                if abs(s) > 0, V = V*(V0/s); end
            end
        end

        % ---- Magnetic impurity 2D: Ns×1 exchange per site -----------
        function VM = mag_2d(p, Nx, Ny)
            Ns = Nx*Ny;
            VM = zeros(Ns,1);
            if isempty(p.imp_pos), return; end
            d   = UQMaT.Disorder.Disorder;
            pos = d.parse_pos_2d(p.imp_pos);
            for k = 1:size(pos,1)
                J  = p.imp_magnetic_val(min(k,numel(p.imp_magnetic_val)));
                VM = VM + d.spread_2d(pos(k,:), J, p.imp_sigma, Nx, Ny, 'gaussian');
            end
        end

        % ---- Parse 2D position input --------------------------------
        function pos = parse_pos_2d(inp)
            inp = inp(:);
            if numel(inp) == 2
                pos = inp(:)';
            elseif mod(numel(inp),2) == 0
                pos = reshape(inp, [], 2);
            else
                error('UQMaT:Disorder:badPos', ...
                    'imp_pos for 2D must be [ix,iy] or Nimp×2 matrix.');
            end
        end

    end   % methods (Static, Access = private)
end
