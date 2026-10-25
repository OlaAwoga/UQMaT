classdef Lanczos < UQMaT.Solver.Base
    % UQMaT.Solver.Diagonalization.Lanczos
    %   Iterative Lanczos (ARPACK). Best for large sparse H, few eigenvalues.
    %
    %   slv = UQMaT.Solver.Diagonalization.Lanczos(sys)
    %   slv = UQMaT.Solver.Diagonalization.Lanczos(sys,'k',10,'sigma',0)
    %
    %   sigma=0 finds eigenvalues nearest to zero (edge/in-gap states).
    %  Citation
    %  --------
    %  Lanczos C., J. Res. Natl. Bur. Stand. 45, 255 (1950). [Lanczos algorithm]
    %  ARPACK (eigs backend): Lehoucq R.B. et al., "ARPACK Users' Guide",
    %    SIAM (1998). [http://www.caam.rice.edu/software/ARPACK]
    properties
        k     (1,1) double {mustBeInteger,mustBePositive} = 20
        sigma  = 'lm'
        tol   (1,1) double = 0
        maxit (1,1) double = 300
    end
    methods
        function obj = Lanczos(sys,varargin)
            obj = obj@UQMaT.Solver.Base(sys);
            for i=1:2:numel(varargin), obj.(varargin{i})=varargin{i+1}; end
        end
        function [E,V] = eigs(obj)
            H   = sparse(obj.sys.H);
            dim = size(H,1);
            k_  = min(obj.k, dim-1);
            opts.tol=obj.tol; opts.maxit=obj.maxit;
            [V,D] = eigs(H, k_, obj.sigma, opts);
            d=real(diag(D)); [E,idx]=sort(d); V=V(:,idx);
        end
        function disp(obj)
            fprintf('  UQMaT.Solver.Diagonalization.Lanczos  [dim=%d k=%d sigma=%s]\n',...
                size(obj.sys.H,1), obj.k, num2str(obj.sigma));
        end
    end
end
