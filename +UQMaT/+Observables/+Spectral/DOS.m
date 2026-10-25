classdef DOS
    % UQMaT.Observables.Spectral.DOS
    %   DOS, LDOS, and spectrum from any UQMaT solver.
    %
    %   obs = UQMaT.Observables.Spectral.DOS(slv)
    %   [E, V]   = obs.spectrum()
    %   [w, D]   = obs.dos()
    %   [w, L]   = obs.ldos()
    %   [w, L]   = obs.ldos(sites)
    %   [w, D_up, D_dn] = obs.dos_spin()    % requires spin=true
    properties (SetAccess=private)
        solver; Ndim (1,1) double
    end
    methods
        function obj = DOS(solver)
            if ~isa(solver,'UQMaT.Solver.Base')
                error('UQMaT:Obs:DOS:badInput','Argument must be UQMaT.Solver.Base.');
            end
            obj.solver = solver;
            obj.Ndim   = solver.sys.Ndim;
        end
        function [E,V] = spectrum(obj), [E,V]=obj.solver.eigs(); end
        function [omega,D] = dos(obj,varargin)
            [omega,eta]=obj.grid(varargin{:});
            E=obj.solver.eigvals(); D=zeros(numel(omega),1);
            for n=1:numel(E), D=D+lor(omega,E(n),eta); end
            D=D/obj.Ndim;
        end
        function [omega,L] = ldos(obj,sites,varargin)
            if nargin<2||isempty(sites), sites=1:obj.Ndim; end
            [omega,eta]=obj.grid(varargin{:});
            [E,V]=obj.solver.eigs();
            Nph=obj.Ndim; Ve=V(1:Nph,:);
            Nw=numel(omega); Ns=numel(sites); L=zeros(Nw,Ns);
            for n=1:numel(E)
                lz=lor(omega,E(n),eta);
                for si=1:Ns
                    L(:,si)=L(:,si)+abs(Ve(sites(si),n))^2*lz;
                end
            end
        end
        function [omega,D_up,D_dn] = dos_spin(obj,varargin)
            if ~obj.solver.sys.params.spin
                error('UQMaT:Obs:DOS:noSpin','spin=true required.'); end
            [omega,eta]=obj.grid(varargin{:});
            [E,V]=obj.solver.eigs(); N=obj.Ndim;
            Nw=numel(omega); D_up=zeros(Nw,1); D_dn=zeros(Nw,1);
            for n=1:numel(E)
                lz=lor(omega,E(n),eta);
                D_up=D_up+sum(abs(V(1:N,n)).^2)*lz;
                D_dn=D_dn+sum(abs(V(N+1:2*N,n)).^2)*lz;
            end
            D_up=D_up/N; D_dn=D_dn/N;
        end
    end
    methods (Access=private)
        function [omega,eta] = grid(obj,varargin)
            p=obj.solver.sys.params; eta=p.eta;
            idx=find(cellfun(@(x)ischar(x)&&strcmpi(x,'omega'),varargin));
            if ~isempty(idx), omega=varargin{idx+1}(:); return; end
            E=obj.solver.sys.E; pad=max(0.2*(E(end)-E(1)),4*eta);
            omega=linspace(E(1)-pad,E(end)+pad,500)';
        end
    end
end
function y=lor(omega,E0,eta), y=(eta/pi)./((omega-E0).^2+eta^2); end
