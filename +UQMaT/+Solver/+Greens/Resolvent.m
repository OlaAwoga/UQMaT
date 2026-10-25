classdef Resolvent < UQMaT.Solver.Base
    % UQMaT.Solver.Greens.Resolvent
    %   G(omega) = ( (omega+i*eta)*I - H )^{-1}  via direct inversion.
    %
    %   slv = UQMaT.Solver.Greens.Resolvent(sys)
    %   [w, D] = slv.dos()
    %   [w, A] = slv.spectral()   % local spectral function (LDOS via G)
    %  Citation
    %  --------
    %  Green's function via resolvent: Economou E.N., "Green's Functions in
    %    Quantum Physics", Springer (2006). [Chapter 2: retarded GF]
    properties
        eta (1,1) double
    end
    methods
        function obj = Resolvent(sys,varargin)
            obj = obj@UQMaT.Solver.Base(sys);
            obj.eta = sys.params.eta;
            for i=1:2:numel(varargin), obj.(varargin{i})=varargin{i+1}; end
        end
        function [E,V] = eigs(obj)
            d = UQMaT.Solver.Diagonalization.Direct(obj.sys);
            [E,V] = d.eigs();
        end
        function G = green(obj, omega_val)
            H = obj.sys.H; dim = size(H,1);
            G = ((omega_val+1i*obj.eta)*eye(dim)-H)\eye(dim);
        end
        function [omega,D] = dos(obj, omega_in)
            omega = obj.grid(omega_in);
            Nph = obj.sys.Ndim; Nw = numel(omega); D = zeros(Nw,1);
            for iw=1:Nw
                G=obj.green(omega(iw));
                D(iw)=-imag(trace(G(1:Nph,1:Nph)))/(pi*Nph);
            end
        end
        function [omega,A] = spectral(obj, omega_in)
            omega=obj.grid(omega_in); Nph=obj.sys.Ndim;
            Nw=numel(omega); A=zeros(Nw,Nph);
            for iw=1:Nw
                G=obj.green(omega(iw));
                for i=1:Nph, A(iw,i)=-imag(G(i,i))/pi; end
            end
        end
        function disp(obj)
            fprintf('  UQMaT.Solver.Greens.Resolvent  [dim=%d eta=%.4g]\n',...
                size(obj.sys.H,1), obj.eta);
        end
    end
    methods (Access=private)
        function omega = grid(obj, omega_in)
            if nargin<2||isempty(omega_in)
                p=obj.sys.params; E=obj.sys.E;
                pad=max(0.2*(E(end)-E(1)),4*obj.eta);
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
