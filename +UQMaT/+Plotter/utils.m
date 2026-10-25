classdef utils
    % UQMaT.Plotter.utils  Shared utilities: colours, axes, option parser, style.
    methods (Static)
        function c=colors(k)
            pal=[0.000,0.447,0.698;0.835,0.369,0.000;0.000,0.620,0.451;
                 0.800,0.475,0.655;0.941,0.894,0.259;0.337,0.706,0.914;0,0,0];
            if nargin==0, c=pal; else, c=pal(mod(k-1,size(pal,1))+1,:); end
        end
        function ax=getAxes(fig,ax_in)
            if ~isempty(ax_in), ax=ax_in; axes(ax);
            elseif ~isempty(fig), figure(fig); ax=gca;
            else, ax=gca; end
        end
        function opt=parseOpts(args,varargin)
            opt=struct();
            for k=1:2:numel(varargin), opt.(varargin{k})=varargin{k+1}; end
            for k=1:2:numel(args)
                if ~isfield(opt,args{k})
                    error('UQMaT:Plotter:unknownOpt','Unknown option: "%s"',args{k}); end
                opt.(args{k})=args{k+1};
            end
        end
        function style(ax)
            set(ax,'FontSize',11,'LineWidth',0.8,'Box','on', ...
                'TickDir','in','XMinorTick','on','YMinorTick','on');
        end

        function c = resolve_colormap(name)
            % Resolve colormap name — handles built-ins and custom names.
            % Custom names not built into base MATLAB R2020b:
            %   'redwhiteblue' | 'rwb' | 'bwr' | 'RdBu'
            %     → built-in red-white-blue (no toolbox needed)
            %   'viridis' | 'magma' | 'inferno' | 'plasma' | 'cividis'
            %     → falls back to 'parula' (MATLAB default)
            switch lower(strtrim(name))
                case {'redwhiteblue','rwb','bwr','rdbu'}
                    n  = 64;
                    r  = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1)];   % blue→white
                    b  = [ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];   % white→red
                    c  = [r; b];
                case {'viridis','magma','inferno','plasma','cividis'}
                    c = parula(128);
                otherwise
                    % Pass through to MATLAB — will error if unknown
                    try
                        c = colormap(name);
                    catch
                        warning('UQMaT:Plotter:badColormap', ...
                            'Unknown colormap "%s" — using parula.', name);
                        c = parula(128);
                    end
            end
        end
    end
end
