classdef (Abstract) Base
    % UQMaT.Models.Base  Abstract base for all UQMaT model classes.
    %
    %   Subclasses must implement buildH() and description().
    %   Models live in:
    %     UQMaT.Models.Chains1D.*      — 1D lattice models
    %     UQMaT.Models.Lattice2D.*     — 2D lattice models
    %     UQMaT.Models.Continuum.*     — continuum models      (planned v3.0)
    %     UQMaT.Models.Quasiclassical.* — quasiclassical models (planned v3.0)

    %  Citation
    %  --------
    %  UQMaT: see CITATION.cff.
    properties (SetAccess = protected)
        params              % UQMaT.Params
        H    (:,:) double   % Hamiltonian matrix
        Ndim (1,1) double   % number of physical sites (electron block)
    end

    methods (Abstract)
        H = buildH(obj)       % construct and return H
        s = description(obj)  % one-line description string
    end

    methods
        function disp(obj)
            fprintf('  %s\n  %s\n  H: %d×%d\n', ...
                class(obj), obj.description(), size(obj.H,1), size(obj.H,2));
        end
    end
end
