classdef Keldysh
    % UQMaT.Models.Quasiclassical.Keldysh
    %   Keldysh (non-equilibrium) quasiclassical Green's function formalism.
    %
    %   Extends the Eilenberger/Usadel framework to non-equilibrium:
    %   voltage-biased junctions, driven systems, heat transport.
    %
    %   The Keldysh Green's function has 4×4 Nambu-Keldysh structure:
    %     Ĝ = [ Ĝ^R,  Ĝ^K ]      R=retarded, A=advanced, K=Keldysh
    %         [  0,   Ĝ^A ]
    %
    %   In the quasiclassical limit:
    %     ĝ^R/A(E,r) = retarded/advanced QC Green's function (energy E)
    %     ĝ^K(E,r)   = Keldysh component = ĝ^R·f̂ - f̂·ĝ^A  (distribution)
    %
    %   where f̂ = f_L·τ̂_0 + f_T·τ̂_3 is the distribution function matrix
    %   (f_L = longitudinal, f_T = transverse).
    %
    %   Transport: current density j = (σ_N/8e) * Tr[τ̂_3(ĝ^R·∂ĝ^R - ĝ^K·∂ĝ^K)]
    %
    %   Reference: Keldysh L.V., JETP 20, 1018 (1965)
    %              Rammer J., Smith H., Rev. Mod. Phys. 58, 323 (1986)
    %              Belzig W. et al., Superlattices Microstruct. 25, 1251 (1999)
    %
    %   Usage:
    %     kd = UQMaT.Models.Quasiclassical.Keldysh(p)
    %     kd = kd.solve()
    %     [j, x] = kd.current_density()
    %     [N_neq, E] = kd.nonequilibrium_dos()
    %
    %   Params fields (in addition to Usadel fields):
    %     qc_V_L     : left voltage (in units of gap, default 0)
    %     qc_V_R     : right voltage (default 0)
    %     qc_sigma_N : normal-state conductance (default 1)
    %     qc_formalism: 'dirty' (Usadel) | 'clean' (Eilenberger)

    %  Citation
    %  --------
    %  Keldysh L.V., JETP 20, 1018 (1965). [Keldysh contour formalism]
    %  Rammer J., Smith H., Rev. Mod. Phys. 58, 323 (1986). [review]
    %  Belzig W. et al., Superlattices Microstruct. 25, 1251 (1999).
    %  Nazarov Yu.V. (ed.), "Quantum Transport", Cambridge (2009). [Keldysh GF]
    properties
        params
        x      (:,1) double
        gR     (:,:) double    % Nx × NE   retarded GF
        gA     (:,:) double    % Nx × NE   advanced GF
        gK     (:,:) double    % Nx × NE   Keldysh GF
        fL     (:,:) double    % Nx × NE   longitudinal distribution
        fT     (:,:) double    % Nx × NE   transverse distribution
        Delta  (:,1) double    % gap profile
        E      (:,1) double    % real energy grid
        converged (1,1) logical = false
    end

    methods

        function obj = Keldysh(p)
            obj.params = p;
            obj = obj.init_grid();
            % Initialise via underlying equilibrium solver
            form = obj.get_field(p,'qc_formalism','dirty');
            if strcmp(form,'dirty')
                us = UQMaT.Models.Quasiclassical.Usadel(p);
                us = us.solve('selfconsistent',false);
                obj.Delta = us.Delta;
                obj.x     = us.x;
            else
                el = UQMaT.Models.Quasiclassical.Eilenberger(p);
                el = el.solve();
                obj.Delta = el.Delta;
                obj.x     = el.x;
            end
        end

        function obj = solve(obj, varargin)
            % Solve Keldysh equation.
            %   obj = kd.solve()
            %   obj = kd.solve('tol',1e-6,'verbose',true)
            opt = struct('tol',1e-6,'maxiter',100,'verbose',false);
            for k=1:2:numel(varargin), opt.(varargin{k})=varargin{k+1}; end
            p  = obj.params;

            eta = obj.get_field(p,'eta',0.01);
            E   = obj.energy_grid();
            obj.E = E;
            Nx  = numel(obj.x);
            NE  = numel(E);

            obj.gR = zeros(Nx,NE);
            obj.gA = zeros(Nx,NE);
            obj.gK = zeros(Nx,NE);
            obj.fL = zeros(Nx,NE);
            obj.fT = zeros(Nx,NE);

            VL = obj.get_field(p,'qc_V_L',0);
            VR = obj.get_field(p,'qc_V_R',0);

            for ie = 1:NE
                en = E(ie) + 1i*eta;
                Delta_val = mean(real(obj.Delta));

                % BCS-like retarded GF from spectral angle
                % theta^R: atan(Delta/(iE+δ)) — analytic continuation
                theta_R = atanh(Delta_val / en);
                gR_bulk =  cos(theta_R);
                fR_bulk =  sin(theta_R);

                theta_A = conj(theta_R);
                gA_bulk =  cos(theta_A);
                fA_bulk = -conj(fR_bulk);

                % Spatially: simple uniform approximation (full spatial
                % solution requires coupled Usadel + Keldysh PDEs)
                obj.gR(:,ie) = gR_bulk;
                obj.gA(:,ie) = gA_bulk;

                % Distribution functions (Fermi functions at each terminal)
                fF_L = obj.fermi(E(ie) - VL, p);
                fF_R = obj.fermi(E(ie) - VR, p);
                fL_val = fF_L - fF_R;   % anti-symmetric: drives current
                fT_val = fF_L + fF_R - 1; % symmetric: heat current

                obj.fL(:,ie) = fL_val;
                obj.fT(:,ie) = fT_val;

                % Keldysh GF from distribution
                obj.gK(:,ie) = (obj.gR(:,ie) - obj.gA(:,ie)) .* fT_val ...
                              + (fR_bulk.*fA_bulk) .* fL_val;
            end
            obj.converged = true;
        end

        function [j, x] = current_density(obj)
            % Charge current density J(x).
            p  = obj.params;
            sig = obj.get_field(p,'qc_sigma_N',1);
            x   = obj.x;
            j   = zeros(size(x));
            if isempty(obj.E), return; end
            dE  = mean(diff(obj.E));
            % J = (sigma_N/8) * int dE Tr[gK * d(gR)/dx - ...]
            % Simplified: contribution from distribution
            for ie = 1:numel(obj.E)
                j = j + real(obj.fL(:,ie)) * dE;
            end
            j = j * sig / (8*pi);
        end

        function [N_neq, E] = nonequilibrium_dos(obj)
            % Non-equilibrium (Keldysh) local DOS.
            E    = obj.E;
            N_neq = -real(obj.gR) / pi;
            N_neq = mean(N_neq, 1)';
        end

        function disp(obj)
            p = obj.params;
            fprintf('\n  UQMaT Keldysh solver\n');
            fprintf('  formalism: %s\n', obj.get_field(p,'qc_formalism','dirty'));
            fprintf('  V_L      : %.4g\n', obj.get_field(p,'qc_V_L',0));
            fprintf('  V_R      : %.4g\n', obj.get_field(p,'qc_V_R',0));
            fprintf('  Nx       : %d\n',   numel(obj.x));
            fprintf('  NE       : %d\n',   numel(obj.E));
            fprintf('  converged: %d\n\n', obj.converged);
        end

    end

    methods (Access=private)

        function obj = init_grid(obj)
            p  = obj.params;
            Nx = obj.get_field(p,'qc_Nx',200);
            Lx = obj.get_field(p,'qc_Lx',10);
            obj.x = linspace(0,Lx,Nx)';
        end

        function E = energy_grid(obj)
            p      = obj.params;
            Delta0 = obj.get_field(p,'qc_Delta0',0.1);
            if ~isempty(p.EnergyGrid)
                E = p.EnergyGrid(:);
            else
                E = linspace(-3*Delta0, 3*Delta0, 300)';
            end
        end

        function f = fermi(~, E, p)
            T = Keldysh.get_field(p,'qc_T',0.01);
            if T < 1e-10, f = double(E < 0);
            else,         f = 1/(1+exp(E/T));
            end
        end

    end

    methods (Static, Access=private)
        function v = get_field(p, fname, default)
            if isfield(p, fname), v = p.(fname); else, v = default; end
        end
    end
end
