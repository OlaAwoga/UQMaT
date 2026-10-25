classdef System
    % UQMaT.System  Routes Params to the correct model and diagonalises.
    %
    %   sys = UQMaT.System(p)      p = UQMaT.Params
    %   sys = UQMaT.System(mdl)    mdl = any UQMaT.Models.Base subclass
    %
    %   Routing:
    %     system='chain'   →  UQMaT.Models.Chains1D.Chain
    %     system='ssh'     →  UQMaT.Models.Chains1D.SSH
    %     system='kitaev'  →  UQMaT.Models.Chains1D.Kitaev
    %     system='square'  →  UQMaT.Models.Lattice2D.Square.Lattice
    %     system='kagome'  →  UQMaT.Models.Lattice2D.Kagome.Lattice
    %
    %   Properties
    %   ----------
    %     H      : Hamiltonian matrix
    %     E      : sorted real eigenvalues
    %     V      : eigenvectors (columns)
    %     params : UQMaT.Params
    %     model  : model object
    %     Ndim   : number of physical sites (electron block)

    %  Citation
    %  --------
    %  UQMaT: see CITATION.cff.
    properties (SetAccess = private)
        H      (:,:) double
        E      (:,1) double
        V      (:,:) double
        params
        model
        Ndim   (1,1) double
    end

    methods
        function obj = System(input)
            if nargin < 1, input = UQMaT.Params(); end

            if isa(input, 'UQMaT.Params')
                obj.params = input;
                obj.model  = UQMaT.System.route(input);
            elseif isa(input, 'UQMaT.Models.Base')
                obj.model  = input;
                obj.params = input.params;
            else
                error('UQMaT:System:badInput', ...
                    'Argument must be UQMaT.Params or UQMaT.Models.Base.');
            end

            obj.H    = obj.model.H;
            obj.Ndim = obj.model.Ndim;

            [V_, E_]    = eig(obj.H, 'vector');
            [obj.E,idx] = sort(real(E_));
            obj.V       = V_(:, idx);
        end

        function disp(obj)
            fprintf('  UQMaT.System\n');
            fprintf('  %s\n', obj.model.description());
            fprintf('  H: %d×%d\n', size(obj.H,1), size(obj.H,2));
            fprintf('  Emin=%.6g  Emax=%.6g\n', obj.E(1), obj.E(end));
            pos = obj.E(obj.E >  1e-10);
            neg = obj.E(obj.E < -1e-10);
            if ~isempty(pos) && ~isempty(neg)
                fprintf('  Gap=%.6g\n', pos(1)-neg(end));
            end
            if obj.params.sc
                UQMaT.SC.BdG.check_ph(obj.H);
            end
        end
    end

    methods (Static, Access = private)
        function mdl = route(p)
            switch p.system
                case 'chain'
                    mdl = UQMaT.Models.Chains1D.Chain(p);
                case 'ssh'
                    mdl = UQMaT.Models.Chains1D.SSH(p);
                case 'kitaev'
                    mdl = UQMaT.Models.Chains1D.Kitaev(p);
                case 'hnkitaev'
                    mdl = UQMaT.Models.Chains1D.HNKitaev(p);
                case 'square'
                    mdl = UQMaT.Models.Lattice2D.Square.Lattice(p);
                case 'kagome'
                    mdl = UQMaT.Models.Lattice2D.Kagome.Lattice(p);
                case 'tmd'
                    mdl = UQMaT.Models.Lattice2D.TMD.Lattice(p);
                case 'cubic'
                    mdl = UQMaT.Models.Lattice3D.Cubic.Lattice(p);
                case 'diamond'
                    mdl = UQMaT.Models.Lattice3D.Diamond.Lattice(p);
                case 'junction1d'
                    mdl = UQMaT.Models.Heterostructure.Junction1D(p);
                case 'chain1d_square2d'
                    mdl = UQMaT.Models.Heterostructure.Chain1D_Square2D(p);
                otherwise
                    error('UQMaT:System:unknown', ...
                        'Unknown system: "%s"', p.system);
            end
        end
    end
end
