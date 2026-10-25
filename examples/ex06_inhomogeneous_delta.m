%% ex06_inhomogeneous_delta.m
%  Inhomogeneous pairing: proximity effect (1D) and vortex core (2D).
%  Each system diagonalised once via UQMaT.EV.

%
%  Solver options (change 'direct' to any of these):
%    ev = UQMaT.EV(sys, 'direct')                        % dense eig()   — default, all states
%    ev = UQMaT.EV(sys, 'lanczos', 'k',20, 'sigma',0)   % sparse eigs() — k states near E=0
%    ev = UQMaT.EV(sys, 'resolvent', 'eta',0.02)         % G(ω) via inversion
%    ev = UQMaT.EV(sys, 'chebyshev', 'M',512)            % KPM — very large systems
uqmat_path();

N   = 80;
t   = 1.0;
mu  = 0.0;
Ds0 = 0.4;

%% ---- 1D proximity effect ----
% Step: SC on right half only
Ds_step        = zeros(N,1);
Ds_step(N/2+1:end) = Ds0;

% Exponential decay from right end
x_sites        = (1:N)';
Ds_decay       = Ds0 * exp(-(N - x_sites) / 15);

% Build SC struct with site-resolved Delta_s
% Build systems and diagonalise with explicit solver
p_uni   = UQMaT.Params('system','chain','N',N,'t',t,'mu',mu,...
    'SC', struct('Delta_s',Ds0,   'sym','onsite'));
p_step  = UQMaT.Params('system','chain','N',N,'t',t,'mu',mu,...
    'SC', struct('Delta_s',Ds_step, 'sym','onsite'));
p_decay = UQMaT.Params('system','chain','N',N,'t',t,'mu',mu,...
    'SC', struct('Delta_s',Ds_decay,'sym','onsite'));

ev_uni   = UQMaT.EV(UQMaT.System(p_uni),   'direct'); % dense eig(), all states
ev_step  = UQMaT.EV(UQMaT.System(p_step),  'direct'); % dense eig(), all states
ev_decay = UQMaT.EV(UQMaT.System(p_decay), 'direct'); % dense eig(), all states

[w, D_uni]   = ev_uni.dos('range',[-3,3]);
[w, D_step]  = ev_step.dos('range',[-3,3]);
[w, D_decay] = ev_decay.dos('range',[-3,3]);
[w, L_step]  = ev_step.ldos([],'range',[-3,3]);
[w, L_decay] = ev_decay.ldos([],'range',[-3,3]);

%% ---- 2D vortex-like gap profile ----
Nx = 20;  Ny = 20;
[IX, IY] = meshgrid(1:Nx, 1:Ny);
cx = (Nx+1)/2;  cy = (Ny+1)/2;
r  = sqrt((IX-cx).^2 + (IY-cy).^2);
Ds_vortex = Ds0 * tanh(r / 3);   % Ny×Nx

p_vx = UQMaT.Params('system','square',...
    'Nx',Nx,'Ny',Ny,'t',t,'mu',mu,...
    'SC', struct('Delta_s',Ds_vortex(:),'sym','onsite'));
ev_vx = UQMaT.EV(UQMaT.System(p_vx), 'direct'); % dense eig(), all states

[w2, D_vx] = ev_vx.dos('range',[-4,4]);
F_vx       = ev_vx.pair_amplitude();
F_vx_2D    = reshape(F_vx, Ny, Nx);

%% ---- Plot ----
figure(1); clf;
set(gcf,'Name','ex06 — inhomogeneous SC');

subplot(2,4,1);
hold on;
plot(x_sites, Ds_step,  'b','LineWidth',2,'DisplayName','Step');
plot(x_sites, Ds_decay, 'r','LineWidth',2,'DisplayName','Decay');
plot(x_sites, repmat(Ds0,N,1),'k:','LineWidth',1,'DisplayName','Uniform');
xlabel('Site'); ylabel('\Delta_s(i)');
title('Pairing profiles'); legend; grid on;

subplot(2,4,2);
UQMaT.Plotter.DOS.compare(w,{D_uni,D_step,D_decay},...
    'labels',{'Uniform','Step','Decay'},'ax',gca,'title','DOS (1D)');

subplot(2,4,3);
UQMaT.Plotter.LDOS.map1D(w,L_step,'ax',gca,'title','LDOS — step proximity');
hold on; xline(N/2,'w--','LineWidth',1.5);

subplot(2,4,4);
UQMaT.Plotter.LDOS.map1D(w,L_decay,'ax',gca,'title','LDOS — decay proximity');

subplot(2,4,5);
imagesc(1:Nx, 1:Ny, Ds_vortex);
axis xy equal tight; colorbar; colormap(gca,parula);
xlabel('x'); ylabel('y');
title('Vortex \Delta_s(r)');

subplot(2,4,6);
UQMaT.Plotter.DOS.compare(w2,{D_vx},'labels',{'Vortex 2D'},'fill',true,...
    'ax',gca,'title','DOS — 2D vortex');

subplot(2,4,7);
imagesc(1:Nx, 1:Ny, abs(F_vx_2D));
axis xy equal tight; colorbar; colormap(gca,hot);
xlabel('x'); ylabel('y');
title('|F(r)| pair amplitude');
