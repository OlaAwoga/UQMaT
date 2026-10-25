classdef Spectral
    % UQMaT.Observables.SC.Spectral
    %   SC-specific observables for systems with sc=true.
    %   Basis: Psi = (c_up, c_dn, c_up+, c_dn+)
    %
    %   obs = UQMaT.Observables.SC.Spectral(slv)
    %   [w, D]           = obs.dos()
    %   [w, D_e, D_h]    = obs.dos_ph()
    %   [w, D_up, D_dn]  = obs.dos_spin()
    %   [w, L]           = obs.ldos(sites)
    %   [w, L_up, L_dn]  = obs.ldos_spin(sites)
    %   [w, L_e, L_h]    = obs.ldos_ph(sites)
    %   [w, Sx,Sy,Sz]    = obs.spin_texture(sites)
    %   F                = obs.pair_amplitude()
    properties (SetAccess=private)
        solver; N (1,1) double
    end
    methods
        function obj = Spectral(solver)
            if ~isa(solver,'UQMaT.Solver.Base')
                error('UQMaT:Obs:SC:badInput','Argument must be UQMaT.Solver.Base.'); end
            if ~solver.sys.params.sc
                error('UQMaT:Obs:SC:notSC','System must have sc=true.'); end
            obj.solver=solver; obj.N=solver.sys.Ndim;
        end
        function [E,V]=spectrum(obj), [E,V]=obj.solver.eigs(); end
        function [omega,D]=dos(obj,varargin)
            [omega,eta]=obj.grid(varargin{:}); E=obj.solver.eigvals();
            D=zeros(numel(omega),1);
            for n=1:numel(E), D=D+lor(omega,E(n),eta); end
            D=D/obj.N;
        end
        function [omega,D_e,D_h]=dos_ph(obj,varargin)
            [omega,eta]=obj.grid(varargin{:}); [E,V]=obj.solver.eigs(); N=obj.N;
            Nw=numel(omega); D_e=zeros(Nw,1); D_h=zeros(Nw,1);
            for n=1:numel(E)
                lz=lor(omega,E(n),eta);
                D_e=D_e+sum(abs(V(1:2*N,n)).^2)*lz;
                D_h=D_h+sum(abs(V(2*N+1:4*N,n)).^2)*lz;
            end
            D_e=D_e/N; D_h=D_h/N;
        end
        function [omega,D_up,D_dn]=dos_spin(obj,varargin)
            [omega,eta]=obj.grid(varargin{:}); [E,V]=obj.solver.eigs(); N=obj.N;
            Nw=numel(omega); D_up=zeros(Nw,1); D_dn=zeros(Nw,1);
            for n=1:numel(E)
                lz=lor(omega,E(n),eta);
                D_up=D_up+sum(abs(V(1:N,n)).^2)*lz;
                D_dn=D_dn+sum(abs(V(N+1:2*N,n)).^2)*lz;
            end
            D_up=D_up/N; D_dn=D_dn/N;
        end
        function [omega,L]=ldos(obj,sites,varargin)
            if nargin<2||isempty(sites), sites=1:obj.N; end
            [omega,eta]=obj.grid(varargin{:}); [E,V]=obj.solver.eigs(); N=obj.N;
            Nw=numel(omega); Ns=numel(sites); L=zeros(Nw,Ns);
            for n=1:numel(E)
                lz=lor(omega,E(n),eta);
                for si=1:Ns
                    i=sites(si);
                    w=abs(V(i,n))^2+abs(V(i+N,n))^2;
                    L(:,si)=L(:,si)+w*lz;
                end
            end
        end
        function [omega,L_up,L_dn]=ldos_spin(obj,sites,varargin)
            if nargin<2||isempty(sites), sites=1:obj.N; end
            [omega,eta]=obj.grid(varargin{:}); [E,V]=obj.solver.eigs(); N=obj.N;
            Nw=numel(omega); Ns=numel(sites); L_up=zeros(Nw,Ns); L_dn=zeros(Nw,Ns);
            for n=1:numel(E)
                lz=lor(omega,E(n),eta);
                for si=1:Ns
                    i=sites(si);
                    L_up(:,si)=L_up(:,si)+abs(V(i,n))^2*lz;
                    L_dn(:,si)=L_dn(:,si)+abs(V(i+N,n))^2*lz;
                end
            end
        end
        function [omega,L_e,L_h]=ldos_ph(obj,sites,varargin)
            if nargin<2||isempty(sites), sites=1:obj.N; end
            [omega,eta]=obj.grid(varargin{:}); [E,V]=obj.solver.eigs(); N=obj.N;
            Nw=numel(omega); Ns=numel(sites); L_e=zeros(Nw,Ns); L_h=zeros(Nw,Ns);
            for n=1:numel(E)
                lz=lor(omega,E(n),eta);
                for si=1:Ns
                    i=sites(si);
                    we=abs(V(i,n))^2+abs(V(i+N,n))^2;
                    wh=abs(V(i+2*N,n))^2+abs(V(i+3*N,n))^2;
                    L_e(:,si)=L_e(:,si)+we*lz;
                    L_h(:,si)=L_h(:,si)+wh*lz;
                end
            end
        end
        function [omega,Sx,Sy,Sz]=spin_texture(obj,sites,varargin)
            if nargin<2||isempty(sites), sites=1:obj.N; end
            [omega,eta]=obj.grid(varargin{:}); [E,V]=obj.solver.eigs(); N=obj.N;
            s=UQMaT.SC.BdG.sigma();
            Nw=numel(omega); Ns=numel(sites);
            Sx=zeros(Nw,Ns); Sy=zeros(Nw,Ns); Sz=zeros(Nw,Ns);
            for n=1:numel(E)
                lz=lor(omega,E(n),eta);
                for si=1:Ns
                    i=sites(si); psi=[V(i,n);V(i+N,n)];
                    Sx(:,si)=Sx(:,si)+real(psi'*s.sx*psi)*lz;
                    Sy(:,si)=Sy(:,si)+real(psi'*s.sy*psi)*lz;
                    Sz(:,si)=Sz(:,si)+real(psi'*s.sz*psi)*lz;
                end
            end
        end
        function F=pair_amplitude(obj)
            % F(i) = sum_{E_n<0} u_n_up(i) * conj(v_n_dn(i))
            [E,V]=obj.solver.eigs(); N=obj.N; F=zeros(N,1);
            for n=find(E<0)'
                u_up=V(1:N,n); v_dn=V(3*N+1:4*N,n);
                F=F+u_up.*conj(v_dn);
            end
        end
    end
    methods (Access=private)
        function [omega,eta]=grid(obj,varargin)
            p=obj.solver.sys.params; eta=p.eta;
            idx=find(cellfun(@(x)ischar(x)&&strcmpi(x,'omega'),varargin));
            if ~isempty(idx), omega=varargin{idx+1}(:); return; end
            E=obj.solver.sys.E; pad=max(0.2*(E(end)-E(1)),4*eta);
            omega=linspace(E(1)-pad,E(end)+pad,500)';
        end
    end
end
function y=lor(omega,E0,eta), y=(eta/pi)./((omega-E0).^2+eta^2); end
