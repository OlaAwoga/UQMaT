classdef LDOS
    % UQMaT.Plotter.LDOS  Plot local density of states.
    methods (Static)
        function ax=map1D(omega,L,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'cmap','hot','clim',[],'title','LDOS','xlabel','Site');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            imagesc(ax,1:size(L,2),omega,L); axis(ax,'xy');
            colorbar(ax); colormap(ax,opt.cmap);
            if ~isempty(opt.clim), caxis(ax,opt.clim); end
            xlabel(ax,opt.xlabel); ylabel(ax,'\omega'); title(ax,opt.title);
            UQMaT.Plotter.utils.style(ax);
        end
        function ax=map2D(Lmap,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'omega',[],'cmap','hot','clim',[],'title','LDOS (spatial)');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            imagesc(ax,abs(Lmap)); axis(ax,'equal','tight');
            colorbar(ax); colormap(ax,opt.cmap);
            if ~isempty(opt.clim), caxis(ax,opt.clim); end
            xlabel(ax,'x'); ylabel(ax,'y');
            if ~isempty(opt.omega), title(ax,sprintf('LDOS  \\omega=%.4g',opt.omega));
            else, title(ax,opt.title); end
            UQMaT.Plotter.utils.style(ax);
        end
        function ax=sites(omega,L,sites,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'labels',{},'lw',1.5,'title','Site-resolved LDOS');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            cols=UQMaT.Plotter.utils.colors(); hold(ax,'on');
            for k=1:min(numel(sites),size(L,2))
                col=cols(mod(k-1,size(cols,1))+1,:);
                lbl=sprintf('Site %d',sites(k));
                if ~isempty(opt.labels)&&k<=numel(opt.labels), lbl=opt.labels{k}; end
                plot(ax,omega,L(:,k),'Color',col,'LineWidth',opt.lw,'DisplayName',lbl);
            end
            xlabel(ax,'\omega'); ylabel(ax,'LDOS'); title(ax,opt.title);
            legend(ax,'Location','best'); grid(ax,'on'); UQMaT.Plotter.utils.style(ax);
        end
        function [ax1,ax2]=spin_map(omega,L_up,L_dn,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'cmap','hot','clim',[]);
            if ~isempty(opt.fig), figure(opt.fig); clf; else, figure; clf; end
            Ns=size(L_up,2);
            subplot(1,2,1); ax1=gca;
            imagesc(ax1,1:Ns,omega,L_up); axis(ax1,'xy'); colorbar(ax1); colormap(ax1,opt.cmap);
            if ~isempty(opt.clim), caxis(ax1,opt.clim); end
            xlabel(ax1,'Site'); ylabel(ax1,'\omega'); title(ax1,'LDOS \uparrow');
            subplot(1,2,2); ax2=gca;
            imagesc(ax2,1:Ns,omega,L_dn); axis(ax2,'xy'); colorbar(ax2); colormap(ax2,opt.cmap);
            if ~isempty(opt.clim), caxis(ax2,opt.clim); end
            xlabel(ax2,'Site'); ylabel(ax2,'\omega'); title(ax2,'LDOS \downarrow');
        end
    end
end
