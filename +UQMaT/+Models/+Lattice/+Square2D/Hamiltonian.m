classdef Hamiltonian < UQMaT.Models.Base
    % UQMaT.Models.Lattice.Square2D.Hamiltonian
    %   2D square lattice tight-binding model.
    %
    %   Site ordering:  n(ix,iy) = (ix-1)*Ny + iy   (column-major)
    %   Total sites:    Ns = Nx * Ny
    %
    %   Supported configurations (all via UQMaT.Params2D):
    %
    %     Normal state                    Ns × Ns
    %     Spinful normal                  2Ns × 2Ns
    %     s-wave BdG (legacy 2-comp)      2Ns × 2Ns
    %     Nambu 4-component BdG           4Ns × 4Ns
    %       Ψ = (c↑, c↓, c↑†, c↓†)ᵀ
    %
    %   Hamiltonian
    %   -----------
    %   H = -t   Σ_{<ij>}  c†_i c_j
    %     - t2   Σ_{<<ij>>} c†_i c_j
    %     - mu   Σ_i        c†_i c_i
    %     + B⃗·σ⃗  (Zeeman)
    %     + i λ_R Σ_{<ij>} c†_i (σ × d_{ij})_z c_j  (Rashba)
    %     + Δ̂ (pairing)
    %
    %   Rashba in 2D:
    %     bond (ix,iy)→(ix+1,iy): H_R += i λ_R σ_y   (x-bond)
    %     bond (ix,iy)→(ix,iy+1): H_R -= i λ_R σ_x   (y-bond)
    %
    %   Usage
    %   -----
    %     p   = UQMaT.Params2D('Nx',20,'Ny',20,'t',1,'mu',-1);
    %     mdl = UQMaT.Models.Lattice.Square2D.Hamiltonian(p);
    %     sys = UQMaT.System(mdl);
    %
    %     % With Nambu BdG:
    %     p   = UQMaT.Params2D('Nx',20,'Ny',20,'pairing','singlet','Delta_s',0.3,'Bz',0.2);
    %     mdl = UQMaT.Models.Lattice.Square2D.Hamiltonian(p);

    methods
        function obj = Hamiltonian(p)
            if ~isa(p,'UQMaT.Params2D')
                error('UQMaT:Square2D:Hamiltonian:badInput', ...
                    'Argument must be a UQMaT.Params2D object.');
            end
            obj.params = p;
            obj.Ndim   = obj.ndim();
            obj.H      = obj.buildH();
        end

        function n = ndim(obj)
            p = obj.params;
            n = p.Nx * p.Ny;
        end

        function s = family(~)
            s = 'lattice';
        end

        function s = description(obj)
            p = obj.params;
            Ns = p.Nx * p.Ny;
            parts = {sprintf('Square2D %d×%d (%d sites) bcx=%s bcy=%s', ...
                p.Nx, p.Ny, Ns, upper(p.bcx), upper(p.bcy))};
            parts{end+1} = sprintf('t=%.3g μ=%.3g', p.t, p.mu);
            if p.t2 ~= 0,  parts{end+1} = sprintf('t2=%.3g', p.t2); end
            if p.nambu
                parts{end+1} = sprintf('Nambu pairing=%s', p.pairing_mode);
            elseif p.bdg
                parts{end+1} = sprintf('BdG Δ=%.3g', mean(p.sc_Delta_s));
            end
            if p.B(3)~=0||p.B(1)~=0||p.B(2)~=0
                parts{end+1} = sprintf('B=[%.2g,%.2g,%.2g]',p.B(1),p.B(2),p.B(3));
            end
            if p.lambda_R ~= 0
                parts{end+1} = sprintf('λ_R=%.3g', p.lambda_R);
            end
            if ~strcmp(p.disorder,'none')
                parts{end+1} = sprintf('disorder=%s W=%.3g', p.disorder, p.W);
            end
            s = strjoin(parts,' | ');
        end

        function H = buildH(obj)
            p  = obj.params;
            Ns = p.Nx * p.Ny;

            % ---- 1. Scalar kinetic + disorder (Ns×Ns) ----
            h0 = obj.kineticMatrix() - p.mu * speye(Ns);
            h0 = h0 + obj.disorderMatrix(Ns);

            % ---- Route by spin/SC mode ----
            if ~p.spin && ~p.bdg && ~p.nambu
                H = full(h0);

            elseif p.spin && ~p.bdg && ~p.nambu
                H = obj.spinfulNormal(h0);

            elseif ~p.spin && p.bdg && ~p.nambu
                % Legacy 2-component s-wave BdG
                Dvec = obj.expandDelta(p.sc_Delta_s, Ns);
                D    = spdiags(Dvec, 0, Ns, Ns);
                H    = full([h0, D; D', -conj(h0)]);

            elseif p.nambu
                H = obj.nambuBdG(h0);

            else
                % spin + sc (non-Nambu)
                H = obj.spinfulBdG(h0);
            end
        end

        function disp(obj)
            p  = obj.params;
            Ns = p.Nx * p.Ny;
            fprintf('\n  UQMaT.Models.Lattice.Square2D.Hamiltonian\n');
            fprintf('  Grid    : %d × %d = %d sites\n', p.Nx, p.Ny, Ns);
            fprintf('  BC      : bcx=%s, bcy=%s\n', upper(p.bcx), upper(p.bcy));
            fprintf('  Normal  : t=%.4g, t2=%.4g, μ=%.4g\n', p.t, p.t2, p.mu);
            if p.nambu
                fprintf('  BdG     : Nambu 4-comp, pairing=%s\n', p.pairing_mode);
                if ismember(p.pairing_mode,{'singlet','mixed'})
                    fprintf('    Δ_s   = %.4g\n', mean(p.sc_Delta_s_s_vec));
                end
                if ismember(p.pairing_mode,{'triplet','mixed'})
                    fprintf('    d-vec = [%.4g, %.4g, %.4g]\n', ...
                        [0,0,0](1), [0,0,0](2), [0,0,0](3));
                end
            elseif p.bdg
                fprintf('  BdG     : 2-comp, Δ=%.4g\n', mean(p.sc_Delta_s));
            end
            if p.spin
                fprintf('  Zeeman  : Bx=%.4g By=%.4g Bz=%.4g\n', p.B(1), p.B(2), p.B(3));
                fprintf('  Rashba  : λ_R=%.4g\n', p.lambda_R);
            end
            if ~strcmp(p.disorder,'none')
                fprintf('  Disorder: %s W=%.4g\n', p.disorder, p.W);
            end
            fprintf('  H size  : %d × %d\n', size(obj.H,1), size(obj.H,2));
            fprintf('\n');
        end
    end

    % ====================================================================
    methods (Access = private)

        % ---- Kinetic matrix (Ns×Ns sparse) ----------------------------
        function T = kineticMatrix(obj)
            p  = obj.params;
            Nx = p.Nx;  Ny = p.Ny;
            Ns = Nx * Ny;
            geo = UQMaT.Models.Lattice.Square2D.Geometry;

            % Use Kronecker structure: T = Tx⊗Iy + Ix⊗Ty
            Tx = UQMaT.Models.Lattice.BoundaryConditions.BC.hop1D(Nx, p.t, 0, p.bcx);
            Ty = UQMaT.Models.Lattice.BoundaryConditions.BC.hop1D(Ny, p.t, 0, p.bcy);
            T  = kron(Tx, speye(Ny)) + kron(speye(Nx), Ty);

            % Next-nearest hopping (diagonal bonds)
            if p.t2 ~= 0
                bonds_nnn = geo.nnn_bonds(Nx, Ny, p.bcx, p.bcy);
                for b = 1:size(bonds_nnn,1)
                    i = bonds_nnn(b,1);  j = bonds_nnn(b,2);
                    T(i,j) = T(i,j) - p.t2;
                    T(j,i) = T(j,i) - p.t2;
                end
            end
        end

        % ---- Rashba hopping matrices (Ns×Ns, x- and y-bonds) ----------
        function [Rx, Ry] = rashbaMatrices(obj)
            % Rx: x-bond Rashba hoppings (Ns×Ns)
            % Ry: y-bond Rashba hoppings (Ns×Ns)
            p  = obj.params;
            Nx = p.Nx;  Ny = p.Ny;
            Ns = Nx * Ny;

            % x-bonds: (ix,iy)→(ix+1,iy), contributes iλ_R σ_y term
            Rx = UQMaT.Models.Lattice.BoundaryConditions.BC.hop1D(Nx, p.lambda_R, 0, p.bcx);
            Rx = kron(Rx, speye(Ny));

            % y-bonds: (ix,iy)→(ix,iy+1), contributes -iλ_R σ_x term
            Ry = UQMaT.Models.Lattice.BoundaryConditions.BC.hop1D(Ny, p.lambda_R, 0, p.bcy);
            Ry = kron(speye(Nx), Ry);
        end

        % ---- Spinful normal-state (2Ns×2Ns) ---------------------------
        function Hs = spinfulNormal(obj, h0)
            p  = obj.params;
            Ns = p.Nx * p.Ny;
            sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
            sz  = sig.sz;  sx = sig.sx;  sy = sig.sy;

            % Spin-diagonal hopping
            Hs = full(kron(eye(2), h0));

            % Zeeman
            if p.B(1)~=0||p.B(2)~=0||p.B(3)~=0
                Hs = Hs + kron(p.B(3)*sz + p.B(1)*sx + p.B(2)*sy, eye(Ns));
            end

            % Rashba
            if p.lambda_R ~= 0
                [Rx, Ry] = obj.rashbaMatrices();
                % x-bond: i λ_R σ_y;  y-bond: -i λ_R σ_x
                Hs = Hs + full( 1i * kron(sy, Rx) - 1i * kron(sx, Ry) );
            end
        end

        % ---- Legacy spinful BdG (spin+sc, non-Nambu) -------------------
        function H = spinfulBdG(obj, h0)
            p  = obj.params;
            Ns = p.Nx * p.Ny;
            sig = UQMaT.Models.Lattice.Nambu.Basis.sigma();
            sy  = sig.sy;

            Hkin = obj.spinfulNormal(h0);
            Dvec = obj.expandDelta(p.sc_Delta_s, Ns);
            D2N  = kron(1i*sy, spdiags(Dvec,0,Ns,Ns));
            H    = full([Hkin, D2N; D2N', -conj(Hkin)]);
        end

        % ---- Full Nambu 4-component BdG (4Ns×4Ns) ----------------------
        function H = nambuBdG(obj, h0)
            p  = obj.params;
            Ns = p.Nx * p.Ny;
            nb  = UQMaT.Models.Lattice.Nambu.Basis;
            sig = nb.sigma();
            sy  = sig.sy;  sx = sig.sx;  sz = sig.sz;

            % Normal block H_N (2Ns×2Ns)
            H_N = full(kron(eye(2), h0));

            % Zeeman
            if p.B(1)~=0||p.B(2)~=0||p.B(3)~=0
                H_N = H_N + kron(p.B(3)*sz + p.B(1)*sx + p.B(2)*sy, eye(Ns));
            end

            % Rashba
            if p.lambda_R ~= 0
                [Rx, Ry] = obj.rashbaMatrices();
                H_N = H_N + full(1i*kron(sy,Rx) - 1i*kron(sx,Ry));
            end

            % Pairing block Δ̂ (2Ns×2Ns)
            Delta_hat = obj.buildPairing2D(Ns, p);

            % Assemble 4Ns×4Ns
            H = nb.assemble(H_N, Delta_hat);
        end

        % ---- Pairing block for 2D (2Ns×2Ns) ---------------------------
        function D = buildPairing2D(obj, Ns, p)
            nb  = UQMaT.Models.Lattice.Nambu.Basis;
            sig = nb.sigma();
            sy  = sig.sy;
            D   = zeros(2*Ns);

            switch p.pairing_mode
                case 'none'

                case 'singlet'
                    Dvec = p.sc_Delta_s_s_vec;   % Ns-vector
                    D    = kron(1i*sy, spdiags(Dvec,0,Ns,Ns));
                    D    = full(D);

                case 'triplet'
                    % p-wave: bond-dependent pairing
                    % d along z: Δ̂_t = dz·σ_z·iσ_y ⊗ T_hop
                    %   = dz * [[0,1];[-1,0]] but split by bond direction
                    [Tx_hop, Ty_hop] = obj.tripletHopMatrices(p);
                    D_trip = obj.tripletPairing2D([0,0,0], Tx_hop, Ty_hop, Ns);
                    D      = D_trip;

                case 'mixed'
                    Dvec   = p.sc_Delta_s_s_vec;
                    D_s    = full(kron(1i*sy, spdiags(Dvec,0,Ns,Ns)));
                    [Tx_hop, Ty_hop] = obj.tripletHopMatrices(p);
                    D_t    = obj.tripletPairing2D([0,0,0], Tx_hop, Ty_hop, Ns);
                    D      = D_s + D_t;
            end
        end

        function [Tx, Ty] = tripletHopMatrices(obj, p)
            Dt = mean(p.sc_Delta_s_t_vec);
            Tx = full(UQMaT.Models.Lattice.BoundaryConditions.BC.hop1D(p.Nx, Dt, 0, p.bcx));
            Tx = kron(Tx, speye(p.Ny));
            Ty = full(UQMaT.Models.Lattice.BoundaryConditions.BC.hop1D(p.Ny, Dt, 0, p.bcy));
            Ty = kron(speye(p.Nx), Ty);
        end

        function D = tripletPairing2D(~, dvec, Tx, Ty, Ns)
            % (d⃗·σ)(iσ_y) decomposed into x and y bond contributions
            %   x-bond: i*sin(kx) → Tx
            %   y-bond: i*sin(ky) → Ty
            % For d = (dx,dy,dz):
            %   d·σ·iσ_y = [ dz,       dx-i*dy ]
            %               [ dx+i*dy,    -dz   ]
            dx = dvec(1);  dy = dvec(2);  dz = dvec(3);
            dsig_isy = [ dz,         dx-1i*dy ;
                         dx+1i*dy,  -dz        ];
            % Combine x and y hoppings (antisymmetric for p-wave)
            T_bond = Tx - Ty;   % relative sign from px±ipy structure
            D = kron(dsig_isy, T_bond);
        end

        % ---- Disorder matrix (Ns×Ns diagonal) -------------------------
        function Vd = disorderMatrix(obj, Ns)
            p = obj.params;
            if p.seed > 0, rng(p.seed); end
            switch p.disorder
                case 'none'
                    Vd = sparse(Ns, Ns);
                case 'anderson'
                    v  = p.W * (rand(Ns,1) - 0.5);
                    Vd = spdiags(v, 0, Ns, Ns);
                case 'binary'
                    v  = p.W * (2*round(rand(Ns,1)) - 1) / 2;
                    Vd = spdiags(v, 0, Ns, Ns);
                case 'impurity'
                    Vd = sparse(Ns, Ns);
                    for k = 1:numel(p.imp_sites)
                        i = p.imp_sites(k);
                        if i >= 1 && i <= Ns
                            Vd(i,i) = p.imp_potential(k);
                        end
                    end
                otherwise
                    Vd = sparse(Ns, Ns);
            end
        end

        % ---- 1D hopping vector (sparse, for Kronecker) ----------------
        function T = hop1D(~, N, t, t2, bc)
            e1 = -t * ones(N-1,1);
            T  = spdiags([e1; 0], -1, N, N) + spdiags([0; e1], 1, N, N);
            if t2 ~= 0
                e2 = -t2 * ones(N-2,1);
                T  = T + spdiags([e2;0;0],-2,N,N) + spdiags([0;0;e2],2,N,N);
            end
            if strcmp(bc,'pbc')
                T(1,N) = -t;  T(N,1) = -t;
                if t2 ~= 0 && N > 3
                    T(1,N-1)=-t2; T(2,N)=-t2;
                    T(N-1,1)=-t2; T(N,2)=-t2;
                end
            end
        end

        % ---- Expand scalar Delta to Ns-vector -------------------------
        function v = expandDelta(~, Delta, Ns)
            if isscalar(Delta)
                v = repmat(Delta, Ns, 1);
            else
                v = Delta(:);
                if numel(v) ~= Ns
                    error('UQMaT:Square2D:badDelta', ...
                        'Delta must be scalar or length Ns=%d.',Ns);
                end
            end
        end
    end
end
