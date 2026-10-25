classdef SpinTexture
    % UQMaT.Plotter.SpinTexture  Plot spin expectation values.
    methods (Static)
        function [ax1,ax2,ax3]=all(omega,Sx,Sy,Sz,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'cmap','redwhiteblue');
            if ~isempty(opt.fig), figure(opt.fig); clf; else, figure; clf; end
            data={Sx,Sy,Sz}; lbls={'<S_x>','<S_y>','<S_z>'}; axs=gobjects(1,3);
            for k=1:3
                subplot(1,3,k); axs(k)=gca; D=data{k};
                cv=max(abs(D(:)));
                imagesc(axs(k),1:size(D,2),omega,D); axis(axs(k),'xy');
                colorbar(axs(k)); colormap(axs(k),opt.cmap); caxis(axs(k),[-cv,cv]);
                xlabel(axs(k),'Site'); ylabel(axs(k),'\omega'); title(axs(k),lbls{k});
                UQMaT.Plotter.utils.style(axs(k));
            end
            ax1=axs(1); ax2=axs(2); ax3=axs(3);
        end
        function ax=component(omega,S,label,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],'cmap','redwhiteblue');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax); cv=max(abs(S(:)));
            imagesc(ax,1:size(S,2),omega,S); axis(ax,'xy');
            colorbar(ax); colormap(ax, UQMaT.Plotter.utils.resolve_colormap(opt.cmap)); caxis(ax,[-cv,cv]);
            xlabel(ax,'Site'); ylabel(ax,'\omega'); title(ax,label);
            UQMaT.Plotter.utils.style(ax);
        end
    end
end
