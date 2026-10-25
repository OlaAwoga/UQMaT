classdef PhaseDiagram
    % UQMaT.Observables.Topological.PhaseDiagram
    %   Scan parameter space and compute gap maps and topological invariants.
    %
    %   pd = UQMaT.Observables.Topological.PhaseDiagram(base_params)
    %   pd = UQMaT.Observables.Topological.PhaseDiagram(base_params,'Npts',40)
    %
    %   ---- Kitaev ------------------------------------------------
    %   [mu_v, Dp_v, nu_map, gap_map] = pd.kitaev_phase(mu_range, Dp_range)
    %     Scans (mu, Delta_p) space; computes winding number and gap.
    %
    %   ---- SSH ---------------------------------------------------
    %   [v_v, w_v, phi_map, gap_map]  = pd.ssh_phase(v_range, w_range)
    %     Scans (v, w) space; computes Zak phase and gap.
    %
    %   ---- Generic gap map ---------------------------------------
    %   [x_v, y_v, gap_map] = pd.gap_map(param1, range1, param2, range2)
    %     Scans any two Params fields and returns the quasiparticle gap.
    %
    %   Example
    %   -------
    %     p0 = UQMaT.Params('system','kitaev','N',40,'t',1,'bc','pbc');
    %     pd = UQMaT.Observables.Topological.PhaseDiagram(p0,'Npts',40);
    %     [mu_v, Dp_v, nu_map, gap_map] = pd.kitaev_phase([-3 3],[0.05 1.5]);

    %  Citation
    %  --------
    %  Phase diagram methods: standard parameter scan; see, e.g.,
    %  Alicea J., Rep. Prog. Phys. 75, 076501 (2012). [Fig. 4, Kitaev phase diagram]
    %  Bernevig B.A., "Topological Insulators and Superconductors", Princeton (2013).
    properties
        Npts (1,1) double {mustBeInteger,mustBePositive} = 40
    end
    properties (SetAccess = private)
        base_params UQMaT.Params
    end

    methods
        function obj = PhaseDiagram(p, varargin)
            if ~isa(p,'UQMaT.Params')
                error('UQMaT:Topo:PD:badInput', ...
                    'First argument must be a UQMaT.Params.');
            end
            obj.base_params = p;
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        function [mu_v, Dp_v, nu_map, gap_map] = kitaev_phase(obj, mu_range, Dp_range)
            mu_v    = linspace(mu_range(1),  mu_range(2),  obj.Npts);
            Dp_v    = linspace(Dp_range(1),  Dp_range(2),  obj.Npts);
            nu_map  = zeros(obj.Npts);
            gap_map = zeros(obj.Npts);
            p0      = obj.base_params;
            fprintf('Kitaev phase diagram (%d×%d)...\n', obj.Npts, obj.Npts);
            for i = 1:obj.Npts
                for j = 1:obj.Npts
                    p   = UQMaT.Params('system','kitaev','N',p0.N,'bc','pbc',...
                              't',p0.t,'mu',mu_v(i),'Delta_p',Dp_v(j));
                    sys = UQMaT.System(p);
                    inv = UQMaT.Observables.Topological.Bulk.Invariants(sys);
                    nu_map(i,j)  = inv.winding();
                    gap_map(i,j) = obj.gap(sys.E);
                end
            end
            fprintf('Done.\n');
        end

        function [v_v, w_v, phi_map, gap_map] = ssh_phase(obj, v_range, w_range)
            v_v     = linspace(v_range(1), v_range(2), obj.Npts);
            w_v     = linspace(w_range(1), w_range(2), obj.Npts);
            phi_map = zeros(obj.Npts);
            gap_map = zeros(obj.Npts);
            p0      = obj.base_params;
            fprintf('SSH phase diagram (%d×%d)...\n', obj.Npts, obj.Npts);
            for i = 1:obj.Npts
                for j = 1:obj.Npts
                    p   = UQMaT.Params('system','ssh','N',p0.N,'bc','pbc',...
                              'v',v_v(i),'w',w_v(j));
                    sys = UQMaT.System(p);
                    inv = UQMaT.Observables.Topological.Bulk.Invariants(sys);
                    phi_map(i,j) = inv.zak();
                    gap_map(i,j) = obj.gap(sys.E);
                end
            end
            fprintf('Done.\n');
        end

        function [x_v, y_v, gap_map] = gap_map(obj, p1, r1, p2, r2)
            x_v     = linspace(r1(1), r1(2), obj.Npts);
            y_v     = linspace(r2(1), r2(2), obj.Npts);
            gap_map = zeros(obj.Npts);
            p0      = obj.base_params;
            fprintf('Gap map (%s vs %s, %d×%d)...\n', p1, p2, obj.Npts, obj.Npts);
            for i = 1:obj.Npts
                for j = 1:obj.Npts
                    p   = UQMaT.Params('system',p0.system,'N',p0.N,'bc',p0.bc,...
                              't',p0.t,'mu',p0.mu,p1,x_v(i),p2,y_v(j));
                    sys = UQMaT.System(p);
                    gap_map(i,j) = obj.gap(sys.E);
                end
            end
            fprintf('Done.\n');
        end
    end

    methods (Access = private)
        function g = gap(~, E)
            pos = E(E >  1e-10);
            neg = E(E < -1e-10);
            if isempty(pos) || isempty(neg), g = 0;
            else, g = pos(1) - neg(end);
            end
        end
    end
end
