classdef Direct < UQMaT.Solver.Base
    % UQMaT.Solver.Diagonalization.Direct
    %   Dense eigensolver — all eigenvalues and eigenvectors.
    %
    %   Uses MATLAB eig() by default, or the compiled MEX C++/LAPACK
    %   solver (uqmat_eig_mex) when available and MEX is enabled.
    %
    %   Usage:
    %     slv = UQMaT.Solver.Diagonalization.Direct(sys)
    %     [E, V] = slv.eigs()
    %
    %   MEX control (global):
    %     UQMaT.Config.use_mex('auto')   % use MEX if compiled (default)
    %     UQMaT.Config.use_mex(true)     % force MEX
    %     UQMaT.Config.use_mex(false)    % force MATLAB eig()
    %
    %  Citation
    %  --------
    %  MATLAB eig: MathWorks, MATLAB Reference (2024).
    %  MEX LAPACK: Anderson E. et al., LAPACK Users' Guide, SIAM (1999).

    methods
        function obj = Direct(sys)
            obj = obj@UQMaT.Solver.Base(sys);
        end

        function [E, V] = eigs(obj)
            H = full(obj.sys.H);

            if UQMaT.Config.should_use_mex('eig')
                % ── MEX C++ LAPACK path ──────────────────────────────
                if UQMaT.Config.is_verbose()
                    fprintf('UQMaT: using MEX eigensolver (uqmat_eig_mex)\n');
                end
                [E, V] = uqmat_eig_mex(H);
                E = E(:);
                % Sort ascending (LAPACK returns sorted, but ensure)
                [E, idx] = sort(real(E));
                V = V(:, idx);
            else
                % ── MATLAB eig fallback ──────────────────────────────
                if UQMaT.Config.is_verbose()
                    fprintf('UQMaT: using MATLAB eig() (MEX not compiled)\n');
                end
                [V, D] = eig(H, 'vector');
                [E, idx] = sort(real(D));
                E = E(:);
                V = V(:, idx);
            end
        end
    end
end
