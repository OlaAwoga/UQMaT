classdef Spectrum
    % UQMaT.Plotter.Spectrum  Plot eigenvalue spectra.
    methods (Static)
        function ax=eigenvalues(sys,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'color',UQMaT.Plotter.utils.colors(1),'title','');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            E=sys.E; plot(ax,1:numel(E),E,'.','Color',opt.color,'MarkerSize',5);
            xlabel(ax,'State index'); ylabel(ax,'Energy'); yline(ax,0,'k--','Alpha',0.4);
            grid(ax,'on');
            if isempty(opt.title), t=sprintf('Spectrum — %s',sys.params.system);
            else, t=opt.title; end
            title(ax,t); UQMaT.Plotter.utils.style(ax);
        end
        function ax=low_energy(sys,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],'nstates',20,...
                'color',UQMaT.Plotter.utils.colors(1),'title','Low-energy spectrum');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            E=sys.E; [~,idx]=sort(abs(E)); E=sort(E(idx(1:min(opt.nstates,numel(E)))));
            bar(ax,1:numel(E),E,'FaceColor',opt.color,'EdgeColor','none');
            yline(ax,0,'k--','Alpha',0.6); xlabel(ax,'Index'); ylabel(ax,'Energy');
            title(ax,opt.title); grid(ax,'on'); UQMaT.Plotter.utils.style(ax);
        end
        function ax=compare(sys_list,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'labels',{},'title','Spectrum comparison');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax); hold(ax,'on');
            cols=UQMaT.Plotter.utils.colors();
            for k=1:numel(sys_list)
                E=sys_list{k}.E; lbl='';
                if ~isempty(opt.labels)&&k<=numel(opt.labels), lbl=opt.labels{k}; end
                plot(ax,1:numel(E),E,'.','Color',cols(mod(k-1,size(cols,1))+1,:),...
                    'MarkerSize',5,'DisplayName',lbl);
            end
            if ~isempty(opt.labels), legend(ax,'Location','best'); end
            yline(ax,0,'k--','Alpha',0.4); xlabel(ax,'State'); ylabel(ax,'Energy');
            title(ax,opt.title); grid(ax,'on'); UQMaT.Plotter.utils.style(ax);
        end
    end
end
