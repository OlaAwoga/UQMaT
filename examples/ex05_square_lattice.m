%% ex05_square_lattice.m
%  2D square lattice: normal state, d-wave SC, band structure, cylinder modes.

%
%  Solver options (change 'direct' to any of these):
%    ev = UQMaT.EV(sys, 'direct')                        % dense eig()   — default, all states
%    ev = UQMaT.EV(sys, 'lanczos', 'k',20, 'sigma',0)   % sparse eigs() — k states near E=0
%    ev = UQMaT.EV(sys, 'resolvent', 'eta',0.02)         % G(ω) via inversion
%    ev = UQMaT.EV(sys, 'chebyshev', 'M',512)            % KPM — very large systems
uqmat_path();

Nx = 16;  Ny = 16;
t  = 1.0; mu = 0.0;

% ---- Normal state — single diagonalisation ----
p_n  = UQMaT.Params('system','square','Nx',Nx,'Ny',Ny,'t',t,'mu',mu,...
                    'bcx','obc','bcy','obc');
sys_n = UQMaT.System(p_n);
ev_n  = UQMaT.EV(sys_n, 'direct');   % dense eig(), all states
[w, D_n] = ev_n.dos('range',[-5,5]);

% ---- d-wave SC — single diagonalisation ----
p_sc  = UQMaT.Params('system','square','Nx',Nx,'Ny',Ny,'t',t,'mu',mu,...
                     'SC',struct('Delta_s',0.3,'sym','dx2y2'),...
                     'bcx','obc','bcy','obc');
sys_sc = UQMaT.System(p_sc);
ev_sc  = UQMaT.EV(sys_sc, 'direct'); % dense eig(), all states

[w, D_sc]  = ev_sc.dos('range',[-5,5]);
[w, L_sc]  = ev_sc.ldos([],'range',[-5,5]);
F          = ev_sc.pair_amplitude();
F2D        = reshape(F, Ny, Nx);

fprintf('Normal H=%dx%d   SC BdG H=%dx%d\n',...
    ev_n.Hdim, ev_n.Hdim, ev_sc.Hdim, ev_sc.Hdim);

% ---- Band structures ----
[kd_kk, Ek_kk, kl, kt] = UQMaT.Models.Lattice2D.Square.BC.bands(...
    UQMaT.System(UQMaT.Params('system','square','Nx',1,'Ny',1,'t',t,'mu',mu,...
    'bcx','pbc','bcy','pbc')), 200);

[kd_cx, Ek_cx, kl_cx, kt_cx] = UQMaT.Models.Lattice2D.Square.BC.bands(...
    UQMaT.System(UQMaT.Params('system','square','Nx',1,'Ny',Ny,'t',t,'mu',mu,...
    'bcx','pbc','bcy','obc')), 200);

% ---- Plot ----
figure(1); clf;
set(gcf,'Name','ex05 — 2D square lattice');

subplot(2,3,1);
UQMaT.Plotter.Bands.plot2D(kd_kk,Ek_kk,kl,kt,'ax',gca,'title','Normal bands Γ→X→M→Γ');

subplot(2,3,2);
UQMaT.Plotter.Bands.plot2D(kd_cx,Ek_cx,kl_cx,kt_cx,'ax',gca,...
    'title','Cylinder (PBC×OBC) modes');

subplot(2,3,3);
UQMaT.Plotter.DOS.compare(w,{D_n,D_sc},...
    'labels',{'Normal','d-wave SC'},'fill',true,'ax',gca,'title','DOS');

subplot(2,3,4);
UQMaT.Plotter.LDOS.map1D(w,L_sc,'ax',gca,'title','LDOS (site-ordered)');

subplot(2,3,5);
UQMaT.Plotter.Lattice.square2D(p_n,'ax',gca,'title',sprintf('%d×%d lattice',Nx,Ny));

subplot(2,3,6);
imagesc(1:Nx, 1:Ny, abs(F2D));
axis xy equal tight; colorbar; colormap(gca,hot);
xlabel('x'); ylabel('y');
title('|F(i)| pair amplitude (d-wave)');
