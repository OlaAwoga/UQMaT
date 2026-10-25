%% ex01_chain_spectrum.m
%  Normal chain vs s-wave superconducting chain.
%  Demonstrates: UQMaT.EV single diagonalisation, DOS, LDOS, spectrum.

%
%  Solver options (change 'direct' to any of these):
%    ev = UQMaT.EV(sys, 'direct')                        % dense eig()   — default, all states
%    ev = UQMaT.EV(sys, 'lanczos', 'k',20, 'sigma',0)   % sparse eigs() — k states near E=0
%    ev = UQMaT.EV(sys, 'resolvent', 'eta',0.02)         % G(ω) via inversion
%    ev = UQMaT.EV(sys, 'chebyshev', 'M',512)            % KPM — very large systems
uqmat_path();

N  = 60;
t  = 1.0;
mu = 0.0;

% ---- Normal chain — diagonalise once ----
p_n  = UQMaT.Params('system','chain','N',N,'t',t,'mu',mu);
sys_n = UQMaT.System(p_n);
ev_n  = UQMaT.EV(sys_n, 'direct');   % dense eig(), all states

[w, D_n] = ev_n.dos('range',[-3,3]);

% ---- s-wave SC chain — diagonalise once ----
p_sc  = UQMaT.Params('system','chain','N',N,'t',t,'mu',mu,...
                     'SC',struct('Delta_s',0.3,'sym','onsite'));
sys_sc = UQMaT.System(p_sc);
ev_sc  = UQMaT.EV(sys_sc, 'direct');  % dense eig(), all states

[w, D_sc]     = ev_sc.dos('range',[-3,3]);
[w, D_e, D_h] = ev_sc.dos_ph('range',[-3,3]);
[w, L_sc]     = ev_sc.ldos([],'range',[-3,3]);
F             = ev_sc.pair_amplitude();

fprintf('Normal:   H = %d×%d\n', ev_n.Hdim, ev_n.Hdim);
fprintf('SC BdG:   H = %d×%d\n', ev_sc.Hdim, ev_sc.Hdim);
pos = ev_sc.E(ev_sc.E > 1e-8);
neg = ev_sc.E(ev_sc.E < -1e-8);
fprintf('SC gap:   %.4f\n', pos(1) - neg(end));

% ---- Plot ----
figure(1); clf;
set(gcf,'Name','ex01 — chain spectrum');

subplot(2,2,1);
UQMaT.Plotter.Spectrum.compare({UQMaT.System(p_n), UQMaT.System(p_sc)}, ...
    'labels',{'Normal','SC BdG'},'ax',gca,'title','Eigenvalue spectra');

subplot(2,2,2);
UQMaT.Plotter.DOS.compare(w,{D_n,D_sc},...
    'labels',{'Normal','SC'},'fill',true,'ax',gca,'title','DOS');

subplot(2,2,3);
UQMaT.Plotter.DOS.particle_hole(w,D_e,D_h,'ax',gca,'title','Particle-hole DOS');

subplot(2,2,4);
UQMaT.Plotter.LDOS.map1D(w,L_sc,'ax',gca,'title','LDOS — SC chain');
