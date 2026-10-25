classdef DOS
    % UQMaT.Plotter.DOS  Plot density of states.
    methods (Static)
        function ax=plot(omega,D,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'color',UQMaT.Plotter.utils.colors(1),'label','','fill',false,...
                'lw',1.5,'title','DOS');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            if opt.fill
                fill(ax,[omega;flipud(omega)],[D;zeros(size(D))],opt.color,...
                    'FaceAlpha',0.25,'EdgeColor','none'); hold(ax,'on');
            end
            plot(ax,omega,D,'Color',opt.color,'LineWidth',opt.lw,'DisplayName',opt.label);
            xlabel(ax,'\omega'); ylabel(ax,'DOS/site'); title(ax,opt.title);
            grid(ax,'on'); if ~isempty(opt.label), legend(ax,'Location','best'); end
            UQMaT.Plotter.utils.style(ax);
        end
        function ax=compare(omega,D_list,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'labels',{},'fill',false,'lw',1.5,'title','DOS comparison');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            cols=UQMaT.Plotter.utils.colors(); hold(ax,'on');
            for k=1:numel(D_list)
                col=cols(mod(k-1,size(cols,1))+1,:); D=D_list{k}; lbl='';
                if ~isempty(opt.labels)&&k<=numel(opt.labels), lbl=opt.labels{k}; end
                if opt.fill, fill(ax,[omega;flipud(omega)],[D;zeros(size(D))],col,...
                        'FaceAlpha',0.2,'EdgeColor','none'); end
                plot(ax,omega,D,'Color',col,'LineWidth',opt.lw,'DisplayName',lbl);
            end
            xlabel(ax,'\omega'); ylabel(ax,'DOS/site'); title(ax,opt.title);
            grid(ax,'on'); if ~isempty(opt.labels), legend(ax,'Location','best'); end
            UQMaT.Plotter.utils.style(ax);
        end
        function ax=spin(omega,D_up,D_dn,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'lw',1.5,'title','Spin-resolved DOS');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            cols=UQMaT.Plotter.utils.colors(); hold(ax,'on');
            plot(ax,omega,D_up,'Color',cols(1,:),'LineWidth',opt.lw,'DisplayName','\uparrow');
            plot(ax,omega,D_dn,'--','Color',cols(2,:),'LineWidth',opt.lw,'DisplayName','\downarrow');
            xlabel(ax,'\omega'); ylabel(ax,'DOS/site'); title(ax,opt.title);
            legend(ax,'Location','best'); grid(ax,'on'); UQMaT.Plotter.utils.style(ax);
        end
        function ax=particle_hole(omega,D_e,D_h,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'lw',1.5,'title','Particle-hole DOS');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            cols=UQMaT.Plotter.utils.colors(); hold(ax,'on');
            plot(ax,omega,D_e,'Color',cols(1,:),'LineWidth',opt.lw,'DisplayName','electron');
            plot(ax,omega,D_h,'--','Color',cols(2,:),'LineWidth',opt.lw,'DisplayName','hole');
            xlabel(ax,'\omega'); ylabel(ax,'DOS/site'); title(ax,opt.title);
            legend(ax,'Location','best'); grid(ax,'on'); UQMaT.Plotter.utils.style(ax);
        end
        function ax=with_std(omega,D,D_std,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'color',UQMaT.Plotter.utils.colors(1),'lw',1.5,'title','Disorder-averaged DOS');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            fill(ax,[omega;flipud(omega)],[D+D_std;flipud(D-D_std)],opt.color,...
                'FaceAlpha',0.25,'EdgeColor','none'); hold(ax,'on');
            plot(ax,omega,D,'Color',opt.color,'LineWidth',opt.lw);
            xlabel(ax,'\omega'); ylabel(ax,'DOS/site'); title(ax,opt.title);
            grid(ax,'on'); UQMaT.Plotter.utils.style(ax);
        end
    end
end

        % ================================================================
        % Self-consistency convergence
        % ================================================================
        function ax = sc_convergence(history, varargin)
            % Plot the self-consistent gap convergence history.
            %   history: struct returned by obs.selfconsistent()
            %
            %   UQMaT.Plotter.DOS.sc_convergence(history)
            opt = UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'title','SC self-consistency convergence','lw',1.5);
            ax = UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            semilogy(ax, 1:numel(history.error), history.error, ...
                'b','LineWidth',opt.lw);
            xlabel(ax,'Iteration'); ylabel(ax,'max|\Delta\Delta|');
            title(ax,opt.title); grid(ax,'on');
            if history.converged
                xline(ax,numel(history.error),'g--','Converged');
            end
            UQMaT.Plotter.utils.style(ax);
        end
