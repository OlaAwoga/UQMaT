classdef Plotter
    % UQMaT.Plotter  Shared utilities for all UQMaT plotting classes.
    %
    %   These are static helpers used internally by:
    %     UQMaT.Plotter.Spectrum
    %     UQMaT.Plotter.DOS
    %     UQMaT.Plotter.LDOS
    %     UQMaT.Plotter.Bands
    %     UQMaT.Plotter.Phase
    %     UQMaT.Plotter.SpinTexture
    %     UQMaT.Plotter.Lattice
    %
    %   They can also be called directly:
    %     cols = UQMaT.Plotter.colors()         % Ncolors × 3 RGB matrix
    %     col  = UQMaT.Plotter.colors(k)        % k-th colour
    %     ax   = UQMaT.Plotter.getAxes(fig, ax) % figure/axes management
    %     opt  = UQMaT.Plotter.parseOpts(args, defaults...)
    %     UQMaT.Plotter.style(ax)               % apply house style

    methods (Static)

        % ================================================================
        % Colour palette
        % ================================================================
        function c = colors(k)
            % Return the UQMaT colour palette as an N×3 RGB matrix,
            % or the k-th colour if k is specified (wraps around).
            %
            % Designed to be:
            %   - Colourblind-friendly (Okabe-Ito palette base)
            %   - Distinguishable in black & white print
            palette = [
                0.000, 0.447, 0.698;   % blue
                0.835, 0.369, 0.000;   % vermillion
                0.000, 0.620, 0.451;   % bluish green
                0.800, 0.475, 0.655;   % reddish purple
                0.941, 0.894, 0.259;   % yellow
                0.337, 0.706, 0.914;   % sky blue
                0.000, 0.000, 0.000;   % black
            ];

            if nargin == 0
                c = palette;
            else
                idx = mod(k-1, size(palette,1)) + 1;
                c   = palette(idx, :);
            end
        end

        % ================================================================
        % Axes management
        % ================================================================
        function ax = getAxes(fig, ax_in)
            % Return axes: use ax_in if given, else create in fig or new figure.
            if ~isempty(ax_in)
                ax = ax_in;
                axes(ax);
            elseif ~isempty(fig)
                figure(fig);
                ax = gca;
            else
                ax = gca;
            end
        end

        % ================================================================
        % Option parser
        % ================================================================
        function opt = parseOpts(args, varargin)
            % Parse Name-Value pairs with defaults.
            %
            %   opt = parseOpts(user_args, 'key1', default1, 'key2', default2, ...)
            %
            % Example:
            %   opt = UQMaT.Plotter.parseOpts(varargin, ...
            %       'fig',   [], ...
            %       'color', 'b', ...
            %       'lw',    1.5);
            %   then use opt.fig, opt.color, opt.lw.

            if mod(numel(varargin),2) ~= 0
                error('UQMaT:Plotter:parseOpts', ...
                    'Defaults must be Name-Value pairs.');
            end

            % Build defaults struct
            opt = struct();
            for k = 1:2:numel(varargin)
                opt.(varargin{k}) = varargin{k+1};
            end

            % Override with user args
            if mod(numel(args),2) ~= 0
                error('UQMaT:Plotter:parseOpts', ...
                    'Arguments must be Name-Value pairs.');
            end
            for k = 1:2:numel(args)
                name = args{k};
                if ~isfield(opt, name)
                    error('UQMaT:Plotter:unknownOpt', ...
                        'Unknown plot option: "%s"', name);
                end
                opt.(name) = args{k+1};
            end
        end

        % ================================================================
        % House style
        % ================================================================
        function style(ax)
            % Apply the UQMaT house style to axes.
            set(ax, ...
                'FontSize',    11, ...
                'LineWidth',   0.8, ...
                'Box',         'on', ...
                'TickDir',     'in', ...
                'XMinorTick',  'on', ...
                'YMinorTick',  'on');
        end

        % ================================================================
        % Quick overview panel for any System
        % ================================================================
        function overview(sys, slv, varargin)
            % One-call overview: spectrum + DOS + LDOS for any 1D system.
            %
            %   UQMaT.Plotter.overview(sys, slv)
            %   UQMaT.Plotter.overview(sys, slv, 'fig', 10)
            opt = UQMaT.Plotter.parseOpts(varargin, 'fig', []);

            if ~isempty(opt.fig), figure(opt.fig); clf;
            else, figure; clf; end

            p = sys.params;

            % Choose the right observable class
            if p.bdg
                obs = UQMaT.Observables.Conventional.Nambu(slv);
            else
                obs = UQMaT.Observables.Conventional.Spectral(slv);
            end

            [w, D]   = obs.dos();
            [w, L]   = obs.ldos();

            subplot(1,3,1);
            UQMaT.Plotter.Spectrum.eigenvalues(sys, 'ax', gca);

            subplot(1,3,2);
            UQMaT.Plotter.DOS.plot(w, D, 'ax', gca, ...
                'title', 'DOS', 'fill', true);

            subplot(1,3,3);
            UQMaT.Plotter.LDOS.map1D(w, L, 'ax', gca, 'title', 'LDOS');

            sgtitle(sprintf('UQMaT overview — %s  N=%d', ...
                p.system, p.N), 'FontSize', 13);
        end

    end
end
