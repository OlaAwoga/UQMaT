classdef Bands
    % UQMaT.Plotter.Bands  Plot band structures.
    %   Uses UQMaT.Models.Lattice2D.Square.BC.bands(sys) for 2D.
    %   Uses UQMaT.Models.Chains1D.BC.bands(sys) for 1D.
    methods (Static)
        function ax=plot1D(kn,Ek,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'title','Band structure','lw',1.5,'gap_line',true,'color',[]);
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            cols=UQMaT.Plotter.utils.colors(); hold(ax,'on');
            for b=1:size(Ek,2)
                if isempty(opt.color), col=cols(mod(b-1,size(cols,1))+1,:);
                else, col=opt.color; end
                plot(ax,kn/pi,Ek(:,b),'Color',col,'LineWidth',opt.lw);
            end
            if opt.gap_line, yline(ax,0,'k--','Alpha',0.4); end
            xlim(ax,[-1,1]); xticks(ax,[-1,-0.5,0,0.5,1]);
            xticklabels(ax,{'-\pi','-\pi/2','0','\pi/2','\pi'});
            xlabel(ax,'k/\pi'); ylabel(ax,'Energy'); title(ax,opt.title);
            grid(ax,'on'); UQMaT.Plotter.utils.style(ax);
        end
        function ax=plot2D(kdist,Ek,klabels,kticks,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'title','Band structure','lw',1.5,'gap_line',true,'color',[]);
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            cols=UQMaT.Plotter.utils.colors(); hold(ax,'on');
            for b=1:size(Ek,2)
                if isempty(opt.color), col=cols(mod(b-1,size(cols,1))+1,:);
                else, col=opt.color; end
                plot(ax,kdist,Ek(:,b),'Color',col,'LineWidth',opt.lw);
            end
            set(ax,'XTick',kdist(kticks),'XTickLabel',klabels);
            for kt=kticks(2:end-1), xline(ax,kdist(kt),'k--','Alpha',0.5); end
            if opt.gap_line, yline(ax,0,'k:','Alpha',0.4); end
            ylabel(ax,'Energy'); title(ax,opt.title);
            grid(ax,'on'); UQMaT.Plotter.utils.style(ax);
        end
        function [ax_b,ax_d]=panel(kn,Ek,omega,D,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'title','Bands & DOS','lw',1.5);
            if ~isempty(opt.fig), figure(opt.fig); clf; else, figure; clf; end
            subplot(1,3,1:2); ax_b=gca;
            UQMaT.Plotter.Bands.plot1D(kn,Ek,'ax',ax_b,'lw',opt.lw,'title',opt.title);
            subplot(1,3,3); ax_d=gca;
            plot(ax_d,D,omega,'b','LineWidth',opt.lw);
            xlabel(ax_d,'DOS'); ylabel(ax_d,'Energy');
            yline(ax_d,0,'k--','Alpha',0.4); grid(ax_d,'on');
            ylim(ax_d,ylim(ax_b)); title(ax_d,'DOS');
            UQMaT.Plotter.utils.style(ax_d);
        end
    end
end
