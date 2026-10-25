classdef Lattice
    % UQMaT.Plotter.Lattice  Visualise lattice geometry.
    methods (Static)
        function ax=chain1D(p,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'highlight',[],'title','');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            N=p.N; x=1:N; hold(ax,'on');
            for i=1:N-1, plot(ax,[x(i),x(i+1)],[0,0],'k-','LineWidth',1.5); end
            if strcmp(p.bc,'pbc')
                ax_=linspace(x(N),x(N)+0.5,5); ay=[0,0.3,0.5,0.3,0];
                plot(ax,[ax_,linspace(x(N)+0.5,(x(1)+x(N))/2,5)],...
                    [ay(1:5),linspace(0.5,0.6,5)],'b--','LineWidth',1.5);
            end
            scatter(ax,x,zeros(1,N),50,'ko','filled');
            if ~isempty(opt.highlight)
                scatter(ax,x(opt.highlight),zeros(1,numel(opt.highlight)),80,'ro','filled');
            end
            ylim(ax,[-0.3,0.8]); axis(ax,'off');
            ttl=opt.title; if isempty(ttl), ttl=sprintf('Chain N=%d bc=%s',N,upper(p.bc)); end
            title(ax,ttl);
        end
        function ax=square2D(p,varargin)
            opt=UQMaT.Plotter.utils.parseOpts(varargin,'fig',[],'ax',[],...
                'highlight',[],'title','');
            ax=UQMaT.Plotter.utils.getAxes(opt.fig,opt.ax);
            Nx=p.Nx; Ny=p.Ny; hold(ax,'on');
            for ix=1:Nx, for iy=1:Ny
                if ix<Nx, plot(ax,[ix,ix+1],[iy,iy],'k-','LineWidth',0.8);
                elseif strcmp(p.bcx,'pbc'), plot(ax,[ix,ix+0.4],[iy,iy],'b--','LineWidth',0.8); end
                if iy<Ny, plot(ax,[ix,ix],[iy,iy+1],'k-','LineWidth',0.8);
                elseif strcmp(p.bcy,'pbc'), plot(ax,[ix,ix],[iy,iy+0.4],'r--','LineWidth',0.8); end
            end, end
            [IX,IY]=meshgrid(1:Nx,1:Ny); scatter(ax,IX(:),IY(:),30,'ko','filled');
            if ~isempty(opt.highlight)
                for n=opt.highlight(:)'
                    [ixh,iyh]=UQMaT.Models.Lattice2D.Square.Geometry.inv_idx(n,Ny);
                    scatter(ax,ixh,iyh,60,'ro','filled');
                end
            end
            axis(ax,'equal','tight'); box(ax,'on');
            xlabel(ax,'x'); ylabel(ax,'y');
            ttl=opt.title;
            if isempty(ttl), ttl=sprintf('%d×%d bcx=%s bcy=%s',Nx,Ny,upper(p.bcx),upper(p.bcy)); end
            title(ax,ttl); UQMaT.Plotter.utils.style(ax);
        end
    end
end
