classdef Invariants
    % UQMaT.Observables.Topological.Invariants
    %   Topological invariants for 1D lattice models.
    %
    %   inv = UQMaT.Observables.Topological.Invariants(sys)
    %   nu  = inv.winding()     % Kitaev: winding number
    %   phi = inv.zak()         % SSH: Zak phase
    %   [M,pf0,pfpi] = inv.pfaffian()  % Kitaev: Majorana number
    properties (SetAccess=private)
        sys UQMaT.System
    end
    properties
        Nk (1,1) double {mustBeInteger,mustBePositive} = 200
    end
    methods
        function obj=Invariants(sys,varargin)
            if ~isa(sys,'UQMaT.System')
                error('UQMaT:Topo:badInput','Argument must be UQMaT.System.'); end
            obj.sys=sys;
            for i=1:2:numel(varargin), obj.(varargin{i})=varargin{i+1}; end
        end
        function phi=zak(obj)
            p=obj.sys.params;
            if ~strcmp(p.system,'ssh'), error('UQMaT:Topo:notSSH','Zak requires system=ssh'); end
            kv=linspace(-pi,pi,obj.Nk+1); uocc=cell(obj.Nk+1,1);
            for ik=1:obj.Nk+1
                k=kv(ik); off=p.v+p.w*exp(-1i*k);
                Hk=[p.mu,off;conj(off),p.mu]; [V,D]=eig(Hk);
                [~,io]=sort(real(diag(D))); uocc{ik}=V(:,io(1));
            end
            M=1;
            for ik=1:obj.Nk, M=M*det(uocc{ik}'*uocc{ik+1}); end
            phi=mod(-imag(log(M))+pi,2*pi)-pi;
        end
        function nu=winding(obj)
            p=obj.sys.params;
            if ~strcmp(p.system,'kitaev'), error('UQMaT:Topo:notKitaev','winding requires kitaev'); end
            kv=linspace(-pi,pi,obj.Nk); kv=kv(1:end-1);
            Dp=mean(p.Delta_p_vec);
            dk=-2*p.t*cos(kv)-p.mu+1i*2*Dp*sin(kv);
            phase=unwrap(angle(dk));
            nu=round((phase(end)-phase(1))/(2*pi));
        end
        function [M_num,pf0,pfpi]=pfaffian(obj)
            p=obj.sys.params;
            if ~strcmp(p.system,'kitaev'), error('UQMaT:Topo:notKitaev','Pfaffian requires kitaev'); end
            H0=obj.kitaevHk(0); Hpi=obj.kitaevHk(pi);
            pf0=real((H0-H0.')/2); pf0=pf0(1,2);
            pfpi=real((Hpi-Hpi.')/2); pfpi=pfpi(1,2);
            M_num=sign(pf0*pfpi);
        end
        function disp(obj)
            p=obj.sys.params;
            fprintf('  UQMaT.Observables.Topological.Invariants  system=%s\n',p.system);
            try
                switch p.system
                    case 'ssh'
                        phi=obj.zak();
                        fprintf('  Zak phase: %.4f rad  (%s)\n',phi,...
                            ternary(abs(phi)>pi/2,'topological','trivial'));
                    case 'kitaev'
                        nu=obj.winding(); [M,~,~]=obj.pfaffian();
                        fprintf('  Winding: %d  (%s)\n',nu,...
                            ternary(nu~=0,'topological','trivial'));
                        fprintf('  Pfaffian: %+d  (%s)\n',M,...
                            ternary(M==-1,'topological','trivial'));
                end
            catch e
                fprintf('  (could not compute: %s)\n',e.message);
            end
        end
    end
    methods (Access=private)
        function Hk=kitaevHk(obj,k)
            p=obj.sys.params; ek=-2*p.t*cos(k)-p.mu; Dp=mean(p.Delta_p_vec);
            Hk=[ek,2i*Dp*sin(k);-2i*Dp*sin(k),-ek];
        end
    end
end
function s=ternary(c,a,b), if c, s=a; else, s=b; end, end
