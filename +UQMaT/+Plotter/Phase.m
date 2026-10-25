classdef Phase
    % UQMaT.Plotter.Phase  Plot topological phase diagrams and gap maps.
    methods (Static)
        function ax=winding_map(x_vec,y_vec,nu_map,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'cmap','redwhiteblue','xlabel','Param 1','ylabel','Param 2','title','Winding number');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            imagesc(ax,y_vec,x_vec,nu_map); axis(ax,'xy');
            colorbar(ax); colormap(ax, UQMaT.Plotter.utils.resolve_colormap(opt.cmap));
            xlabel(ax,opt.xlabel); ylabel(ax,opt.ylabel); title(ax,opt.title);
            UQMaT.Plotter.utils.style(ax);
        end
        function ax=gap_map(x_vec,y_vec,gap,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'cmap','parula','xlabel','Param 1','ylabel','Param 2','title','Gap');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            imagesc(ax,y_vec,x_vec,gap); axis(ax,'xy');
            colorbar(ax); colormap(ax, UQMaT.Plotter.utils.resolve_colormap(opt.cmap)); hold(ax,'on');
            contour(ax,y_vec,x_vec,gap,[0.01 0.01],'w--','LineWidth',1.5);
            xlabel(ax,opt.xlabel); ylabel(ax,opt.ylabel); title(ax,opt.title);
            UQMaT.Plotter.utils.style(ax);
        end
        function [ax1,ax2]=combined(x_vec,y_vec,inv_map,gap,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],...
                'xlabel','Param 1','ylabel','Param 2');
            if ~isempty(opt.fig), figure(opt.fig); clf; else, figure; clf; end
            subplot(1,2,1); ax1=gca;
            UQMaT.Plotter.Phase.winding_map(x_vec,y_vec,inv_map,'ax',ax1,...
                'xlabel',opt.xlabel,'ylabel',opt.ylabel,'title','Invariant');
            subplot(1,2,2); ax2=gca;
            UQMaT.Plotter.Phase.gap_map(x_vec,y_vec,gap,'ax',ax2,...
                'xlabel',opt.xlabel,'ylabel',opt.ylabel,'title','Gap');
        end
        function ax=ipr(E,ipr_vals,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'cmap','hot','title','IPR');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            scatter(ax,E,ipr_vals,8,ipr_vals,'filled');
            colorbar(ax); colormap(ax, UQMaT.Plotter.utils.resolve_colormap(opt.cmap));
            xlabel(ax,'Energy'); ylabel(ax,'IPR'); title(ax,opt.title);
            grid(ax,'on'); UQMaT.Plotter.utils.style(ax);
        end
    end
end
