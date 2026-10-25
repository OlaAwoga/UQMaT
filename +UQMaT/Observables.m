classdef Observables
    % UQMaT.Observables  Legacy convenience wrapper (v1.0/v1.1 compatibility).
    %
    %   Wraps UQMaT.Observables.Conventional.Spectral using the Direct solver.
    %   New code should use the explicit solver + observable classes:
    %
    %     slv = UQMaT.Solver.Diagonalization.Direct(sys)
    %     obs = UQMaT.Observables.Conventional.Spectral(slv)
    %
    %   This class remains for backward compatibility and quick interactive use.

    properties (SetAccess = private)
        inner   % UQMaT.Observables.Conventional.Spectral
        sys     UQMaT.System
    end

    methods
        function obj = Observables(sys)
            obj.sys   = sys;
            slv       = UQMaT.Solver.Diagonalization.Direct(sys);
            obj.inner = UQMaT.Observables.Conventional.Spectral(slv);
        end

        function [E, V]   = spectrum(obj),          [E,V]   = obj.inner.spectrum(); end
        function [w, D]   = dos(obj, varargin),     [w,D]   = obj.inner.dos(varargin{:}); end
        function [w, L]   = ldos(obj, varargin),    [w,L]   = obj.inner.ldos(varargin{:}); end

        function disp(obj)
            fprintf('  UQMaT.Observables (legacy wrapper → Conventional.Spectral)\n');
        end
    end
end
