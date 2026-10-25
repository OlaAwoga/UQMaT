classdef Lattice < UQMaT.Models.Base
    % UQMaT.Models.Lattice2D.Square.Lattice
    %   2D square tight-binding lattice.
    %   Site index: n(ix,iy) = (ix-1)*Ny + iy  (column-major)
    %   Total sites: Ns = Nx * Ny
    %
    %   H sizes:
    %     sc=false, spin=false  →  Ns × Ns
    %     sc=false, spin=true   →  2Ns × 2Ns
    %     sc=true               →  4Ns × 4Ns  Ψ=(c↑,c↓,c↑†,c↓†)ᵀ
    %
    %   Independent BC per direction:
    %     bcx='obc', bcy='obc'  →  finite slab
    %     bcx='pbc', bcy='obc'  →  cylinder (periodic in x)
    %     bcx='obc', bcy='pbc'  →  cylinder (periodic in y)
    %     bcx='pbc', bcy='pbc'  →  torus (fully periodic)

    %  Citation
    %  --------
    %  Tight-binding square lattice: Ashcroft N.W., Mermin N.D., "Solid State
    %    Physics", Holt (1976). [Chapter 10: square lattice dispersion]
    %  Tinkham M., "Introduction to Superconductivity", McGraw-Hill (1996).
    %    [Chapter 7: BdG on lattice, d-wave pairing]
    methods
        function obj = Lattice(p)
            obj.params = p;
            obj.Ndim   = p.Nx * p.Ny;
            obj.H      = obj.buildH();
        end

        function s = description(obj)
            p  = obj.params;
            Ns = p.Nx*p.Ny;
            if p.sc
                s = sprintf('Square2D SC BdG 4×%d (%d×%d) bcx=%s bcy=%s', ...
                    Ns, p.Nx, p.Ny, upper(p.bcx), upper(p.bcy));
            elseif p.spin
                s = sprintf('Square2D spinful 2×%d (%d×%d) bcx=%s bcy=%s', ...
                    Ns, p.Nx, p.Ny, upper(p.bcx), upper(p.bcy));
            else
                s = sprintf('Square2D normal %d (%d×%d) bcx=%s bcy=%s', ...
                    Ns, p.Nx, p.Ny, upper(p.bcx), upper(p.bcy));
            end
        end

        function H = buildH(obj)
            p   = obj.params;
            Nx  = p.Nx; Ny = p.Ny; Ns = Nx*Ny;
            bdg = UQMaT.SC.BdG;
            hop = @UQMaT.Models.Chains1D.Chain.hop1D;

            % Validate sc
            if p.sc && all(p.sc_Delta_s == 0) && p.sc_Delta_t == 0
                error('UQMaT:Square:noAmplitude', ...
                    ['sc=true but all pairing amplitudes are zero.\n' ...
                     'Set Delta_s, dvec, or Delta_t to a nonzero value.\n' ...
                     'Example: UQMaT.Params(''system'',''square'',''Delta_s'',0.3)']);
            end

            % Order parameter contributions (all channels)
            [D_s_ops, D_t_ops, h_AM, h_CDW] = UQMaT.OrderParameterSymmetry.OPS.build(p);

            % Scalar kinetic + CDW (Kronecker structure)
            Tx = hop(Nx, p.tx, p.t2, p.bcx);
            Ty = hop(Ny, p.ty, p.t2, p.bcy);
            h0 = full(kron(Tx, speye(Ny)) + kron(speye(Nx), Ty)) ...
                 - p.mu*eye(Ns) + p.onsite_val*eye(Ns) ...
                 + diag(UQMaT.Disorder.Disorder.impurity_2d_or_bulk(p, Ns, Nx, Ny)) ...
                 + h_CDW;

            if p.sc
                s   = bdg.sigma();
                H_N = kron(eye(2), h0);
                if any(p.B ~= 0)
                    H_N = H_N + bdg.zeeman(p.B,Ns);
                end
                if p.lambda_R~=0
                    Ax  = full(kron(UQMaT.SC.BdG.rashba_hop(Nx,p.bcx), speye(Ny)));
                    Ay  = full(kron(speye(Nx), UQMaT.SC.BdG.rashba_hop(Ny,p.bcy)));
                    H_N = H_N + p.lambda_R*(1i*kron(s.sy,Ax) - 1i*kron(s.sx,Ay));
                end
                H_N = H_N + h_AM;
                D   = D_s_ops + D_t_ops;
                H   = bdg.assemble(H_N, D);

            elseif p.spin
                s  = bdg.sigma();
                H  = kron(eye(2), h0);
                if any(p.B ~= 0)
                    H = H + bdg.zeeman(p.B,Ns);
                end
                if p.lambda_R~=0
                    Ax = full(kron(UQMaT.SC.BdG.rashba_hop(Nx,p.bcx), speye(Ny)));
                    Ay = full(kron(speye(Nx), UQMaT.SC.BdG.rashba_hop(Ny,p.bcy)));
                    H  = H + p.lambda_R*(1i*kron(s.sy,Ax) - 1i*kron(s.sx,Ay));
                end
                H = H + h_AM;

            else
                % Normal spinless: Ns × Ns
                H = h0;
            end
        end
    end

end
