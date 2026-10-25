classdef (Abstract) Base
    % UQMaT.Solver.Base  Abstract base for all solvers.
    %  Citation
    %  --------
    %  UQMaT: see CITATION.cff.
    properties (SetAccess = protected)
        sys  UQMaT.System
    end
    methods
        function obj = Base(sys)
            if ~isa(sys,'UQMaT.System')
                error('UQMaT:Solver:badInput','Argument must be UQMaT.System.');
            end
            obj.sys = sys;
        end
        function E = eigvals(obj)
            [E,~] = obj.eigs();
        end
    end
    methods (Abstract)
        [E, V] = eigs(obj)
    end
end
