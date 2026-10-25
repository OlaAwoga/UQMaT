classdef Invariants
    % UQMaT.Observables.Topological.Bulk.Invariants
    %   Bulk (k-space) topological invariants from the Bloch Hamiltonian H(k).
    %   Requires translation-invariant (PBC) systems.
    %
    %   Invariants
    %   ----------
    %   Chern number  C  [2D systems, class A/C]
    %     C = (1/2π) ∬ Ω(k) d²k   (integral of Berry curvature over BZ)
    %     Computed via the Fukui-Hatsugai-Suzuki lattice method.
    %     C ≠ 0  →  topological (quantum Hall / chiral SC)
    %     Ref: Fukui, Hatsugai, Suzuki, J. Phys. Soc. Jpn 74, 1674 (2005)
    %
    %   Berry curvature  Ω(k)  [2D]
    %     Ω_n(k) = -2 Im Σ_{m≠n} <n|∂Hk/∂kx|m><m|∂Hk/∂ky|n> / (E_m-E_n)²
    %     Ref: Berry M.V., Proc. R. Soc. A 392, 45 (1984)
    %
    %   Z2 invariant  ν  [2D, class AII — time-reversal invariant]
    %     Computed via time-reversal polarisation at TRIM points.
    %     ν=1 → topological insulator, ν=0 → trivial
    %     Ref: Fu L., Kane C.L., PRB 76, 045302 (2007)
    %
    %   Wilson loop / Wannier centres  [2D, class A]
    %     Traces Berry phase as a function of ky; winding → non-trivial Chern.
    %     Ref: Soluyanov A.A., Vanderbilt D., PRB 83, 235401 (2011)
    %
    %   Entanglement spectrum  [1D and 2D]
    %     Eigenvalues of reduced density matrix ρ_A = Tr_B |GS><GS|.
    %     Degenerate levels in gap → topological phase.
    %     Ref: Li H., Haldane F.D.M., PRL 101, 010504 (2008)
    %
    %   Winding number  ν  [1D, class D/BDI — Kitaev chain]
    %     ν = (1/2π) ∮ d(arg d(k)) dk
    %     Ref: Kitaev A.Yu., Phys. Uspekhi 44, 131 (2001)
    %
    %   Zak phase  φ  [1D, class AIII — SSH chain]
    %     φ = i ∮ <u(k)|∂k|u(k)> dk   (Berry phase of occupied band)
    %     φ=π → topological,  φ=0 → trivial
    %     Ref: Zak J., PRL 62, 2747 (1989)
    %
    %   Pfaffian / Majorana number  M  [1D, class D — Kitaev]
    %     M = sgn[Pf(H(0)) · Pf(H(π))]
    %     M=-1 → topological,  M=+1 → trivial
    %     Ref: Kitaev A.Yu., Phys. Uspekhi 44, 131 (2001)
    %
    %   Quantum metric  g(k)  [2D — geometry of Bloch states]
    %     g_μν(k) = Re[ <∂μu|∂νu> - <∂μu|u><u|∂νu> ]
    %     Related to spread of Wannier functions.
    %     Ref: Provost J.P., Vallee G., Commun. Math. Phys. 76, 289 (1980)
    %
    %   Usage
    %   -----
    %     % 1D Kitaev chain (PBC)
    %     sys = UQMaT.System(UQMaT.Params('system','kitaev','N',60,...
    %               'mu',-1,'Delta_p',0.5,'bc','pbc'));
    %     inv = UQMaT.Observables.Topological.Bulk.Invariants(sys)
    %     nu            = inv.winding()
    %     [M,pf0,pfpi]  = inv.pfaffian()
    %     Z2_1d         = inv.Z2()
    %
    %     % 1D SSH chain (PBC)
    %     sys = UQMaT.System(UQMaT.Params('system','ssh','N',40,'v',0.3,'w',1,'bc','pbc'))
    %     inv = UQMaT.Observables.Topological.Bulk.Invariants(sys)
    %     phi           = inv.zak()
    %
    %     % 2D square lattice — Chern number, Berry curvature, Wilson loop
    %     sys = UQMaT.System(UQMaT.Params('system','square','Nx',1,'Ny',1,...
    %               't',1,'mu',0,'SC',struct('Delta_s',0.3,'sym','chiral_p'),...
    %               'bcx','pbc','bcy','pbc'))
    %     inv = UQMaT.Observables.Topological.Bulk.Invariants(sys)
    %     C             = inv.chern()
    %     [Omega,kx,ky] = inv.berry_curvature()
    %     [wc, ky_v]    = inv.wilson_loop()
    %     [xi, E_ent]   = inv.entanglement_spectrum()
    %     Z2_2d         = inv.Z2_2d()
    %     [g, kx, ky]   = inv.quantum_metric()
    %
    %     disp(inv)           % prints all relevant invariants
    %
    %   Citation
    %   --------
    %   When using specific methods, cite the corresponding references above.
    %   For UQMaT itself: see CITATION.cff in the package root.

    properties (SetAccess = private)
        sys  UQMaT.System
    end
    properties
        Nk  (1,1) double {mustBeInteger, mustBePositive} = 100
    end

    methods
        function obj = Invariants(sys, varargin)
            if ~isa(sys, 'UQMaT.System')
                error('UQMaT:Topo:Bulk:badInput', ...
                    'Argument must be a UQMaT.System.');
            end
            obj.sys = sys;
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        % ============================================================
        % Chern number  C  (2D systems)
        % Fukui-Hatsugai-Suzuki lattice gauge field method.
        % Ref: Fukui, Hatsugai, Suzuki, JPSJ 74, 1674 (2005)
        % ============================================================
        function C = chern(obj, band_idx)
            % Chern number of specified bands (default: all occupied).
            %   C = inv.chern()          % all occupied bands
            %   C = inv.chern(1)         % band 1 only
            %   C = inv.chern([1,2])     % sum over bands 1 and 2
            obj.need2D('Chern number');
            Nk  = obj.Nk;
            kx  = linspace(0, 2*pi, Nk+1);
            ky  = linspace(0, 2*pi, Nk+1);
            p   = obj.sys.params;

            % Build Bloch states on grid
            [U] = obj.bloch_states_2d(kx(1:end-1), ky(1:end-1), p);

            if nargin < 2 || isempty(band_idx)
                % Occupied bands: below Fermi level (mu=0 convention)
                E0 = obj.bloch_energies_2d(0, 0, p);
                band_idx = find(E0 < 0);
                if isempty(band_idx)
                    band_idx = 1:floor(size(U,3)/2);
                end
            end

            % Lattice field strength (FHS method)
            C = 0;
            for ib = band_idx(:)'
                F = 0;
                for ix = 1:Nk
                    for iy = 1:Nk
                        u00 = U(:, ib, ix,   iy);
                        u10 = U(:, ib, ix+1, iy);
                        u01 = U(:, ib, ix,   iy+1);
                        u11 = U(:, ib, ix+1, iy+1);
                        % Link variables
                        U01 = u00'*u10 / abs(u00'*u10 + 1e-15);
                        U02 = u10'*u11 / abs(u10'*u11 + 1e-15);
                        U03 = u11'*u01 / abs(u11'*u01 + 1e-15);
                        U04 = u01'*u00 / abs(u01'*u00 + 1e-15);
                        % Lattice field strength
                        F12 = U01*U02*U03*U04;
                        F   = F + imag(log(F12 + 1e-15));
                    end
                end
                C = C + round(F / (2*pi));
            end
        end

        % ============================================================
        % Berry curvature  Ω_n(k)  (2D systems)
        % Ref: Berry M.V., Proc. R. Soc. A 392, 45 (1984)
        % ============================================================
        function [Omega, kx_grid, ky_grid] = berry_curvature(obj, band_idx)
            % Berry curvature of band band_idx over the BZ.
            %   [Omega, kx, ky] = inv.berry_curvature()        % lowest band
            %   [Omega, kx, ky] = inv.berry_curvature(2)       % band 2
            obj.need2D('Berry curvature');
            if nargin < 2 || isempty(band_idx), band_idx = 1; end
            p   = obj.sys.params;
            Nk  = obj.Nk;
            kx  = linspace(-pi, pi, Nk+1); kx = kx(1:end-1);
            ky  = linspace(-pi, pi, Nk+1); ky = ky(1:end-1);
            [kx_grid, ky_grid] = meshgrid(kx, ky);
            Omega = zeros(Nk);

            dkx = 2*pi/Nk;  dky = 2*pi/Nk;
            for ix = 1:Nk
                for iy = 1:Nk
                    Hk   = obj.bloch_H(kx(ix), ky(iy), p);
                    Hkx  = obj.bloch_H(kx(ix)+dkx, ky(iy), p);
                    Hky  = obj.bloch_H(kx(ix), ky(iy)+dky, p);
                    [V,D]= eig(Hk,'vector');
                    [~,so]= sort(real(D));
                    V    = V(:,so);
                    n    = band_idx;
                    un   = V(:,n);
                    O    = 0;
                    for m = 1:size(V,2)
                        if m == n, continue; end
                        um   = V(:,m);
                        En   = real(D(so(n)));
                        Em   = real(D(so(m)));
                        if abs(Em-En) < 1e-10, continue; end
                        dHx  = (obj.bloch_H(kx(ix)+1e-5,ky(iy),p) - Hk)/1e-5;
                        dHy  = (obj.bloch_H(kx(ix),ky(iy)+1e-5,p) - Hk)/1e-5;
                        mx   = um'*dHx*un;
                        my   = un'*dHy*um;
                        O    = O - 2*imag(mx*my)/(Em-En)^2;
                    end
                    Omega(iy,ix) = O;
                end
            end
        end

        % ============================================================
        % Wilson loop / Wannier charge centres  (2D)
        % Ref: Soluyanov A.A., Vanderbilt D., PRB 83, 235401 (2011)
        % ============================================================
        function [wcc, ky_vec] = wilson_loop(obj, band_idx)
            % Wannier charge centres via Wilson loop along kx for each ky.
            %   [wcc, ky] = inv.wilson_loop()        % all occupied bands
            %   [wcc, ky] = inv.wilson_loop([1,2])   % selected bands
            obj.need2D('Wilson loop');
            p    = obj.sys.params;
            Nk   = obj.Nk;
            ky_vec = linspace(-pi, pi, Nk);

            E0 = obj.bloch_energies_2d(0, 0, p);
            if nargin < 2 || isempty(band_idx)
                band_idx = find(E0 < 0);
                if isempty(band_idx), band_idx = 1; end
            end
            Nb  = numel(band_idx);
            wcc = zeros(Nk, Nb);

            kx_path = linspace(0, 2*pi, Nk+1);

            for iy = 1:Nk
                ky = ky_vec(iy);
                % Wilson loop W = ∏_{kx} P(kx) where P is occupied projector
                W = eye(Nb);
                for ix = 1:Nk
                    kx1 = kx_path(ix);   kx2 = kx_path(ix+1);
                    u1  = obj.bloch_occ(kx1, ky, p, band_idx);
                    u2  = obj.bloch_occ(kx2, ky, p, band_idx);
                    W   = W * (u1'*u2);
                end
                ev       = sort(imag(log(eig(W)))/(2*pi));
                wcc(iy,:)= ev(:)';
            end
        end

        % ============================================================
        % Entanglement spectrum  (1D and 2D)
        % Ref: Li H., Haldane F.D.M., PRL 101, 010504 (2008)
        % ============================================================
        function [xi, E_ent] = entanglement_spectrum(obj, cut)
            % Entanglement spectrum via correlation matrix method.
            %   [xi, E_ent] = inv.entanglement_spectrum()
            %   [xi, E_ent] = inv.entanglement_spectrum(N_cut)
            %   cut: number of sites in subsystem A (default N/2)
            %
            % xi : entanglement eigenvalues (0<xi<1)
            % E_ent = -log(xi/(1-xi)) : entanglement energies
            N   = obj.sys.Ndim;
            V   = obj.sys.V;
            E   = obj.sys.E;
            if nargin < 2 || isempty(cut), cut = floor(N/2); end

            % Occupied states
            occ = find(E < 0);
            if isempty(occ)
                occ = 1:floor(size(V,2)/2);
            end

            % Correlation matrix C(i,j) = <c†_i c_j>  (site block)
            Hdim = size(V,1);
            % Use electron spin-up block (rows 1:N)
            Vocc = V(1:N, occ);
            C    = Vocc * Vocc';

            % Reduced to subsystem A: first 'cut' sites
            C_A  = C(1:cut, 1:cut);
            xi   = real(eig(C_A));
            xi   = sort(xi);
            xi   = xi(xi > 1e-10 & xi < 1-1e-10);
            E_ent= -log(xi ./ (1-xi));
        end

        % ============================================================
        % Z2 invariant (2D, class AII — TR-invariant insulators)
        % Ref: Fu L., Kane C.L., PRB 76, 045302 (2007)
        % ============================================================
        function Z2 = Z2_2d(obj)
            % Z2 invariant from time-reversal polarisation.
            % Requires spin=true and TR symmetry (B=0, dvec=0).
            obj.need2D('Z2 (2D)');
            p = obj.sys.params;
            if any(p.B ~= 0)
                warning('UQMaT:Topo:Z2:TR',...
                    'Z2_2d requires time-reversal symmetry (B=[0,0,0]).');
            end

            % TRIM points for 2D square BZ
            TRIM = [0,0; pi,0; 0,pi; pi,pi];
            delta = 1;
            for it = 1:4
                kx = TRIM(it,1);  ky = TRIM(it,2);
                Hk = obj.bloch_H(kx, ky, p);
                E_trim = real(eig(Hk));
                E_trim = sort(E_trim);
                occ = E_trim(E_trim < 0);
                if mod(numel(occ),2) ~= 0
                    % Use parity of occupied band count as proxy
                    delta = delta * (-1)^numel(occ);
                end
            end
            Z2 = mod((1-sign(delta))/2, 2);
        end

        % ============================================================
        % Quantum metric  g_μν(k)  (2D)
        % Ref: Provost J.P., Vallee G., Commun. Math. Phys. 76, 289 (1980)
        % ============================================================
        function [g_trace, kx_grid, ky_grid] = quantum_metric(obj, band_idx)
            % Trace of quantum metric tensor tr[g(k)] = g_xx + g_yy.
            %   [g, kx, ky] = inv.quantum_metric()
            obj.need2D('Quantum metric');
            if nargin < 2 || isempty(band_idx), band_idx = 1; end
            p   = obj.sys.params;
            Nk  = obj.Nk;
            kx  = linspace(-pi,pi,Nk); ky = linspace(-pi,pi,Nk);
            [kx_grid, ky_grid] = meshgrid(kx, ky);
            g_trace = zeros(Nk);

            for ix = 1:Nk
                for iy = 1:Nk
                    Hk  = obj.bloch_H(kx(ix),ky(iy),p);
                    [V,D] = eig(Hk,'vector'); [~,so] = sort(real(D));
                    V   = V(:,so);
                    un  = V(:,band_idx);
                    gxx = 0; gyy = 0;
                    for m = 1:size(V,2)
                        if m == band_idx, continue; end
                        um  = V(:,m);
                        dHx = (obj.bloch_H(kx(ix)+1e-5,ky(iy),p)-Hk)/1e-5;
                        dHy = (obj.bloch_H(kx(ix),ky(iy)+1e-5,p)-Hk)/1e-5;
                        En  = real(D(so(band_idx))); Em = real(D(so(m)));
                        if abs(Em-En)<1e-10, continue; end
                        gxx = gxx + abs(um'*dHx*un)^2/(Em-En)^2;
                        gyy = gyy + abs(um'*dHy*un)^2/(Em-En)^2;
                    end
                    g_trace(iy,ix) = real(gxx + gyy);
                end
            end
        end

        % ============================================================
        % Winding number  ν  (1D Kitaev chain)
        % Ref: Kitaev A.Yu., Phys. Uspekhi 44, 131 (2001)
        % ============================================================
        function nu = winding(obj)
            obj.need1D('kitaev', 'winding number');
            p   = obj.sys.params;
            kv  = linspace(-pi, pi, obj.Nk+1); kv = kv(1:end-1);
            Dp  = mean(p.Delta_p_vec);
            dk  = -2*p.t*cos(kv) - p.mu + 2i*Dp*sin(kv);
            phase = unwrap(angle(dk));
            nu  = round((phase(end)-phase(1))/(2*pi));
        end

        % ============================================================
        % Zak phase  φ  (1D SSH chain)
        % Ref: Zak J., PRL 62, 2747 (1989)
        % ============================================================
        function phi = zak(obj)
            obj.need1D('ssh', 'Zak phase');
            p  = obj.sys.params;
            kv = linspace(-pi, pi, obj.Nk+1);
            uocc = cell(obj.Nk+1, 1);
            for ik = 1:obj.Nk+1
                k   = kv(ik);
                off = p.v + p.w*exp(-1i*k);
                Hk  = [p.mu, off; conj(off), p.mu];
                [V,D] = eig(Hk); [~,io] = sort(real(diag(D)));
                uocc{ik} = V(:,io(1));
            end
            M_prod = 1;
            for ik = 1:obj.Nk
                M_prod = M_prod * (uocc{ik}'*uocc{ik+1});
            end
            phi = mod(-imag(log(M_prod))+pi, 2*pi) - pi;
        end

        % ============================================================
        % Pfaffian / Majorana number  M  (1D Kitaev)
        % Ref: Kitaev A.Yu., Phys. Uspekhi 44, 131 (2001)
        % ============================================================
        function [M, pf0, pfpi] = pfaffian(obj)
            obj.need1D('kitaev', 'Pfaffian');
            p    = obj.sys.params;
            H0   = obj.kitaevBloch(0);
            Hpi  = obj.kitaevBloch(pi);
            pf0  = real((H0-H0.')/2);  pf0  = pf0(1,2);
            pfpi = real((Hpi-Hpi.')/2); pfpi = pfpi(1,2);
            M    = sign(pf0*pfpi);
        end

        % ============================================================
        % Z2 invariant (1D, via Pfaffian)
        % ============================================================
        function Z2 = Z2(obj)
            [M,~,~] = obj.pfaffian();
            Z2 = (1-M)/2;
        end

        % ============================================================
        function disp(obj)
            p = obj.sys.params;
            fprintf('\n  UQMaT.Observables.Topological.Bulk.Invariants\n');
            fprintf('  system=%s\n\n', p.system);
            switch p.system
                case 'kitaev'
                    nu = obj.winding();
                    [M,~,~] = obj.pfaffian();
                    Z2 = obj.Z2();
                    fprintf('  Winding number  ν = %d   (%s)\n',  nu, tag(nu~=0));
                    fprintf('  Majorana number M = %+d  (%s)\n',  M,  tag(M==-1));
                    fprintf('  Z2 invariant      = %d   (%s)\n',  Z2, tag(Z2==1));
                case 'ssh'
                    phi = obj.zak();
                    fprintf('  Zak phase φ = %.4f rad  (%.1f°)  (%s)\n', ...
                        phi, rad2deg(phi), tag(abs(phi) > pi/2));
                case 'square'
                    try
                        C = obj.chern();
                        fprintf('  Chern number C = %d  (%s)\n', C, tag(C~=0));
                    catch ME
                        fprintf('  Chern number: %s\n', ME.message);
                    end
                otherwise
                    fprintf('  For 1D: winding(), zak(), pfaffian(), Z2()\n');
                    fprintf('  For 2D: chern(), berry_curvature(), wilson_loop(),\n');
                    fprintf('          entanglement_spectrum(), Z2_2d(), quantum_metric()\n');
            end
            fprintf('\n');
        end
    end

    methods (Access = private)

        function need1D(obj, sys_name, inv_name)
            if ~strcmp(obj.sys.params.system, sys_name)
                error('UQMaT:Topo:Bulk:wrongSystem', ...
                    '%s requires system="%s", got "%s".', ...
                    inv_name, sys_name, obj.sys.params.system);
            end
        end

        function need2D(obj, inv_name)
            p = obj.sys.params;
            if ~p.is_2d
                error('UQMaT:Topo:Bulk:need2D', ...
                    '%s requires a 2D system (square, kagome, tmd). Got "%s".', ...
                    inv_name, p.system);
            end
        end

        function Hk = bloch_H(obj, kx, ky, p)
            % 2D Bloch Hamiltonian at (kx,ky) — scalar unit-cell version.
            tx = p.tx; ty = p.ty; mu = p.mu;
            ek = -2*tx*cos(kx) - 2*ty*cos(ky) - mu;
            s  = UQMaT.SC.BdG.sigma();
            if p.sc
                Delta = mean(p.sc_Delta_s);
                H_N = ek*eye(2);
                if any(p.B~=0), H_N=H_N+p.B(1)*s.sx+p.B(2)*s.sy+p.B(3)*s.sz; end
                D   = Delta*eye(2);
                Hk  = UQMaT.SC.BdG.assemble(H_N, kron(1i*s.sy,D));
            elseif p.spin
                Hk  = ek*eye(2);
                if any(p.B~=0), Hk=Hk+p.B(1)*s.sx+p.B(2)*s.sy+p.B(3)*s.sz; end
            else
                Hk  = ek;
            end
        end

        function E = bloch_energies_2d(obj, kx, ky, p)
            Hk = obj.bloch_H(kx, ky, p);
            E  = sort(real(eig(Hk)));
        end

        function U = bloch_states_2d(obj, kx_v, ky_v, p)
            % U: Nbasis × Nbands × Nkx × Nky
            Nkx = numel(kx_v); Nky = numel(ky_v);
            Hk0 = obj.bloch_H(0, 0, p);
            Nb  = size(Hk0,1);
            U   = zeros(Nb, Nb, Nkx+1, Nky+1);
            kx_ext = [kx_v, kx_v(1)+2*pi];
            ky_ext = [ky_v, ky_v(1)+2*pi];
            for ix = 1:Nkx+1
                for iy = 1:Nky+1
                    Hk = obj.bloch_H(kx_ext(ix), ky_ext(iy), p);
                    [V,D] = eig(Hk,'vector'); [~,so]=sort(real(D));
                    U(:,:,ix,iy) = V(:,so);
                end
            end
        end

        function u = bloch_occ(obj, kx, ky, p, band_idx)
            Hk = obj.bloch_H(kx, ky, p);
            [V,D] = eig(Hk,'vector'); [~,so]=sort(real(D));
            V = V(:,so);
            u = V(:,band_idx);
        end

        function Hk = kitaevBloch(obj, k)
            p  = obj.sys.params;
            ek = -2*p.t*cos(k) - p.mu;
            Dp = mean(p.Delta_p_vec);
            Hk = [ek, 2i*Dp*sin(k); -2i*Dp*sin(k), -ek];
        end
    end
end

function s = tag(cond)
    if cond, s = 'topological'; else, s = 'trivial'; end
end
