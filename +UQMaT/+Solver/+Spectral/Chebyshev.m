classdef Chebyshev < UQMaT.Solver.Base
    % UQMaT.Solver.Spectral.Chebyshev
    %   Kernel Polynomial Method (KPM) — DOS and LDOS without diagonalisation.
    %   O(N*M) cost. Best for large N where direct eig is too slow.
    %
    %   slv = UQMaT.Solver.Spectral.Chebyshev(sys)
    %   slv = UQMaT.Solver.Spectral.Chebyshev(sys,'M',512)
    %   [w,D] = slv.dos()
    %   [w,L] = slv.ldos(sites)
    %  Citation
    %  --------
    %  KPM (kernel polynomial method): Silver R.N., Roder H., PRE 56, 4822 (1997).
    %  Weisse A. et al., Rev. Mod. Phys. 78, 275 (2006). [comprehensive review]
    properties
        M      (1,1) double {mustBeInteger,mustBePositive} = 256
        Nrand  (1,1) double {mustBeInteger,mustBeNonnegative} = 10
        kernel (1,:) char   = 'jackson'
    end
    properties (SetAccess=private)
        Hsc (:,:) double; a (1,1) double; b (1,1) double
    end
    methods
        function obj = Chebyshev(sys,varargin)
            obj = obj@UQMaT.Solver.Base(sys);
            for i=1:2:numel(varargin), obj.(varargin{i})=varargin{i+1}; end
            obj = obj.rescale();
        end
        function [E,V] = eigs(obj)
            d=UQMaT.Solver.Diagonalization.Direct(obj.sys); [E,V]=d.eigs();
        end
        function [omega,D] = dos(obj,omega_in)
            omega=obj.grid(omega_in); dim=size(obj.Hsc,1);
            Nph=obj.sys.Ndim; g=obj.kernelW(); mu=zeros(obj.M,1);
            if obj.Nrand>0
                for ir=1:obj.Nrand
                    r=rand(dim,1)*2-1; r=r/norm(r); mu=mu+obj.moments(r);
                end
                mu=mu/obj.Nrand;
            else
                for i=1:dim; e=zeros(dim,1); e(i)=1; mu=mu+obj.moments(e); end
                mu=mu/dim;
            end
            D=obj.reconstruct(omega,mu.*g)/Nph;
        end
        function [omega,L] = ldos(obj,sites,omega_in)
            if nargin<2||isempty(sites), sites=1:obj.sys.Ndim; end
            if nargin<3, omega_in=[]; end
            omega=obj.grid(omega_in); g=obj.kernelW();
            dim=size(obj.Hsc,1); Nw=numel(omega); Ns=numel(sites);
            L=zeros(Nw,Ns);
            for si=1:Ns
                e=zeros(dim,1); e(sites(si))=1;
                mu=obj.moments(e);
                L(:,si)=obj.reconstruct(omega,mu.*g);
            end
        end
        function disp(obj)
            fprintf('  UQMaT.Solver.Spectral.Chebyshev  [dim=%d M=%d kernel=%s]\n',...
                size(obj.sys.H,1),obj.M,obj.kernel);
        end
    end
    methods (Access=private)
        function obj=rescale(obj)
            E=obj.sys.E; Emax=E(end)*1.05; Emin=E(1)*1.05;
            obj.a=(Emax-Emin)/2; obj.b=(Emax+Emin)/2;
            obj.Hsc=(obj.sys.H-obj.b*speye(size(obj.sys.H,1)))/obj.a;
        end
        function mu=moments(obj,v)
            M=obj.M; mu=zeros(M,1); a0=v; a1=obj.Hsc*v;
            mu(1)=v'*a0; if M>1, mu(2)=v'*a1; end
            for m=3:M; a2=2*obj.Hsc*a1-a0; mu(m)=v'*a2; a0=a1; a1=a2; end
        end
        function g=kernelW(obj)
            M=obj.M; n=(0:M-1)';
            if strcmp(obj.kernel,'jackson')
                g=((M-n+1).*cos(pi*n/(M+1))+sin(pi*n/(M+1)).*cot(pi/(M+1)))/(M+1);
            else
                lam=4; g=sinh(lam*(1-n/M))/sinh(lam);
            end
        end
        function D=reconstruct(obj,omega,gmu)
            E_sc=max(-1+1e-8,min(1-1e-8,(omega-obj.b)/obj.a));
            D=zeros(size(omega));
            for m=1:obj.M, D=D+gmu(m)*cos((m-1)*acos(E_sc)); end
            D=D*2./(pi*obj.a.*sqrt(1-E_sc.^2));
        end
        function omega=grid(obj,omega_in)
            if nargin<2||isempty(omega_in)
                p=obj.sys.params; E=obj.sys.E;
                pad=max(0.1*(E(end)-E(1)),0.1);
                if ~isempty(p.EnergyGrid)
                omega=p.EnergyGrid(:);
            else
                omega=linspace(E(1)-pad,E(end)+pad,500)';
            end
            else
                omega=omega_in(:);
            end
        end
    end
end
