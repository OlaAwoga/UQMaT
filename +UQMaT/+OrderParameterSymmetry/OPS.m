classdef OPS
    % UQMaT.OrderParameterSymmetry.OPS
    %   Builds all order parameter contributions from Params structs.
    %
    %   Interface (set in Params):
    %     SC  = struct('Delta_s',0.3, 'sym','onsite',   % singlet
    %                  'Delta_t',0.1, 'sym_t','px')     % triplet
    %
    %     AM  = struct('J',0.5,'sym','gwave')            % single component
    %     AM  = struct('J_s',0.3,'sym_s','onsite', ...   % multiple components
    %                  'J_dxy',0.2,'sym_dxy','dxy')
    %     AM  = {struct('J',0.3,'sym','onsite'), ...     % cell array
    %            struct('J',0.2,'sym','dxy')}
    %
    %     CDW = struct('Delta',0.2,'sym','onsite')
    %
    %   Usage in model buildH:
    %     [D_s, D_t, h_AM, h_CDW] = UQMaT.OrderParameterSymmetry.OPS.build(p)
    %
    %   Print all available symmetries:
    %     UQMaT.OrderParameterSymmetry.OPS.describe()

    %  Citation
    %  --------
    %  Singlet/triplet classification: Sigrist M., Ueda K., Rev. Mod. Phys. 63,
    %    239 (1991). [gap symmetry in unconventional SC]
    %  Altermagnetism: Smejkal L. et al., Science Adv. 8, eabl4539 (2022).
    %  CDW: Gruner G., Rev. Mod. Phys. 60, 1129 (1988).
    methods (Static)

        function [D_s, D_t, h_AM, h_CDW] = build(p)
            % Returns:
            %   D_s    : 2N×2N singlet pairing block
            %   D_t    : 2N×2N triplet pairing block
            %   h_AM   : 2N×2N altermagnet spin block (sum of all components)
            %   h_CDW  : N×N   CDW site-diagonal block

            % Heterostructures handle OPS internally — return zeros
            if ismember(p.system,{'junction1d','chain1d_square2d'})
                D_s=0; D_t=0; h_AM=0; h_CDW=0; return;
            end

            bdg = UQMaT.SC.BdG;
            s   = bdg.sigma();
            N   = UQMaT.OrderParameterSymmetry.OPS.Nsites(p);

            % ---- Singlet SC ----
            if any(p.sc_Delta_s ~= 0)
                sym = UQMaT.OrderParameterSymmetry.OPS.get_sym(p, p.sc_sym_s);
                if all(all(sym.D == diag(diag(sym.D))))
                    % On-site (diagonal): use site-resolved Delta_s_vec
                    D_s = kron(1i*s.sy, diag(p.Delta_s_vec(1:N)));
                else
                    % Bond-resolved: uniform amplitude × bond matrix
                    D_s = kron(1i*s.sy, diag(p.Delta_s_vec(1:size(sym.D,1))) * sym.D);
                end
            else
                D_s = zeros(2*N);
            end

            % ---- Triplet SC ----
            if p.sc_Delta_t ~= 0
                sym = UQMaT.OrderParameterSymmetry.OPS.get_sym(p, p.sc_sym_t);
                D_t = kron(1i*s.sy, p.sc_Delta_t * sym.D);
            else
                D_t = zeros(2*N);
            end

            % ---- Altermagnet (sum of all components) ----
            h_AM = zeros(2*N);
            for k = 1:numel(p.am_components)
                c   = p.am_components{k};
                sym = UQMaT.OrderParameterSymmetry.OPS.get_sym(p, c.sym);
                % Each component: J * σ_z ⊗ diag(sublattice pattern of sym)
                h_AM = h_AM + c.J * kron(s.sz, diag(sym.sublattice));
            end

            % ---- CDW ----
            if p.cdw_Delta ~= 0
                sym   = UQMaT.OrderParameterSymmetry.OPS.get_sym(p, p.cdw_sym);
                sw    = sym.site_weights;
                if all(sw == 0), sw = sym.sublattice; end
                h_CDW = p.cdw_Delta * diag(sw);
            else
                h_CDW = zeros(N);
            end
        end

        function sym = get_sym(p, label)
            label = lower(strtrim(label));
            switch p.system
                case {'chain','kitaev'}
                    sym = UQMaT.OrderParameterSymmetry.Chain1D.get(label, p.N, p.bc);
                case 'ssh'
                    sym = UQMaT.OrderParameterSymmetry.Chain1D.get(label, 2*p.N, p.bc);
                case 'square'
                    sym = UQMaT.OrderParameterSymmetry.Square2D.get(label, p.Nx, p.Ny, p.bcx, p.bcy);
                case 'honeycomb'
                    sym = UQMaT.OrderParameterSymmetry.Honeycomb.get(label, p.Nx, p.Ny, p.bc);
                case 'kagome'
                    sym = UQMaT.OrderParameterSymmetry.Kagome.get(label, p.Nx, p.Ny, p.bc);
                case 'tmd'
                    sym = UQMaT.OrderParameterSymmetry.Square2D.get(label, p.Nx, p.Ny, p.bcx, p.bcy);
                case {'cubic','diamond'}
                    % 3D: use onsite symmetry as default
                    Ns = p.Nx*p.Ny*p.Nz;
                    if strcmp(p.system,'diamond'), Ns = 2*Ns; end
                    sym = struct('D',eye(Ns),'site_weights',ones(Ns,1),...
                                 'sublattice',(-1).^(1:Ns)','label',label,'system',p.system);
                case {'junction1d','chain1d_square2d'}
                    % Heterostructures: return zero OPS (handled internally)
                    N = 1;
                    sym = struct('D',0,'site_weights',0,'sublattice',0,'label','none','system','hetero');
                otherwise
                    error('UQMaT:OPS:unknownSystem', ...
                        'No OPS class for system "%s"', p.system);
            end
        end

        function list = available(system)
            switch lower(system)
                case {'chain','ssh','kitaev'}
                    list = UQMaT.OrderParameterSymmetry.Chain1D.available();
                case 'square'
                    list = UQMaT.OrderParameterSymmetry.Square2D.available();
                case 'honeycomb'
                    list = UQMaT.OrderParameterSymmetry.Honeycomb.available();
                case 'kagome'
                    list = UQMaT.OrderParameterSymmetry.Kagome.available();
                otherwise
                    list = {};
            end
        end

        function describe()
            fprintf('\n  UQMaT Order Parameter Symmetry\n');
            fprintf('  ==============================\n\n');
            fprintf('  SC  struct fields:\n');
            fprintf('    Delta_s  : singlet amplitude\n');
            fprintf('    sym      : singlet symmetry (default: ''onsite'')\n');
            fprintf('    Delta_t  : triplet amplitude\n');
            fprintf('    sym_t    : triplet symmetry (default: ''px'')\n\n');
            fprintf('  AM  struct fields (multiple components supported):\n');
            fprintf('    J    + sym              : single component\n');
            fprintf('    J_s  + sym_s, J_dxy + sym_dxy, ...  : multi-component\n');
            fprintf('    {struct(''J'',..,''sym'',..), ...}      : cell array\n\n');
            fprintf('  CDW struct fields:\n');
            fprintf('    Delta + sym\n\n');
            fprintf('  Available symmetries per lattice:\n');
            for sys = {'chain','square','honeycomb'}
                list = UQMaT.OrderParameterSymmetry.OPS.available(sys{1});
                fprintf('    %-10s: %s\n', sys{1}, strjoin(list,', '));
            end
            fprintf('\n  Examples:\n');
            fprintf("    p = UQMaT.Params('system','chain','N',60, ...\n");
            fprintf("        'SC', struct('Delta_s',0.2,'sym','onsite', ...\n");
            fprintf("                    'Delta_t',0.1,'sym_t','px'), ...\n");
            fprintf("        'AM', struct('J',0.5,'sym','gwave'))\n\n");
            fprintf("    p = UQMaT.Params('system','square','Nx',20,'Ny',20, ...\n");
            fprintf("        'SC',  struct('Delta_s',0.3,'sym','dx2y2'), ...\n");
            fprintf("        'AM',  struct('J_s',0.3,'sym_s','onsite','J_dxy',0.2,'sym_dxy','dxy'), ...\n");
            fprintf("        'CDW', struct('Delta',0.1,'sym','onsite'))\n\n");
        end

    end

    methods (Static, Access = private)
        function N = Nsites(p)
            if p.is_2d
                N = p.Nx * p.Ny;
            elseif strcmp(p.system,'ssh')
                N = 2*p.N;
            else
                N = p.N;
            end
        end
    end
end
