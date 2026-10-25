%% ex02_topological_sc.m
%  Topological SC: Rashba chain + Zeeman + s-wave pairing.
%  Demonstrates: MZM edge states, spin-resolved LDOS, Majorana polarization.

%
%  Solver options (change 'direct' to any of these):
%    ev = UQMaT.EV(sys, 'direct')                        % dense eig()   — default, all states
%    ev = UQMaT.EV(sys, 'lanczos', 'k',20, 'sigma',0)   % sparse eigs() — k states near E=0
%    ev = UQMaT.EV(sys, 'resolvent', 'eta',0.02)         % G(ω) via inversion
%    ev = UQMaT.EV(sys, 'chebyshev', 'M',512)            % KPM — very large systems
uqmat_path();

N        = 80;
t        = 1.0;
mu       = 0.0;
Delta_s  = 0.3;
lambda_R = 0.5;
B_topo   = [0, 0, 0.6];   % Bz > Delta_s → topological
B_triv   = [0, 0, 0.1];   % Bz < Delta_s → trivial

% ---- Topological phase — single diagonalisation ----
p_t  = UQMaT.Params('system','chain','N',N,'t',t,'mu',mu,...
                    'SC',struct('Delta_s',Delta_s,'sym','onsite'),...
                    'B',B_topo,'lambda_R',lambda_R);
sys_t = UQMaT.System(p_t);
ev_t  = UQMaT.EV(sys_t, 'direct');   % dense eig(), all states

[w, D_t]        = ev_t.dos('range',[-2,2],'Npts',400);
[w, L_t]        = ev_t.ldos([],'range',[-2,2],'Npts',400);
[w, D_up, D_dn] = ev_t.dos_spin('range',[-2,2],'Npts',400);
F_t             = ev_t.pair_amplitude();

% In-gap states
idx_t = ev_t.edge_idx();
fprintf('Topological: %d in-gap states  E = %s\n',...
    numel(idx_t), num2str(ev_t.E(idx_t)','%.5f '));

% Majorana polarization
mp_t = ev_t.mp();
[MP_L, MP_R, MP_LR] = mp_t.left_right(idx_t);
fprintf('MP_L=%.4f  MP_R=%.4f  MP_LR=%.4f\n', MP_L, MP_R, MP_LR);

% ---- Trivial phase ----
p_v  = UQMaT.Params('system','chain','N',N,'t',t,'mu',mu,...
                    'SC',struct('Delta_s',Delta_s,'sym','onsite'),...
                    'B',B_triv,'lambda_R',lambda_R);
sys_v = UQMaT.System(p_v);
ev_v  = UQMaT.EV(sys_v, 'direct');   % dense eig(), all states
[w, D_v] = ev_v.dos('range',[-2,2],'Npts',400);

% ---- Band structure ----
sys_k = UQMaT.System(UQMaT.Params('system','chain','N',1,'t',t,'mu',mu,...
    'SC',struct('Delta_s',Delta_s,'sym','onsite'),'B',B_topo,'lambda_R',lambda_R));
[kn, Ek] = UQMaT.Models.Chains1D.ChainBC.bands(sys_k, 400);

% ---- Real-space invariants (from cached EV) ----
inv_t    = ev_t.realspace_invariants();
[rho, x] = inv_t.edge_density();

% ---- Plot ----
figure(1); clf;
set(gcf,'Name','ex02 — topological SC');

subplot(2,3,1);
UQMaT.Plotter.DOS.compare(w,{D_t,D_v},...
    'labels',{'Topo (Bz>Δ)','Trivial (Bz<Δ)'},'ax',gca,...
    'title','DOS: topo vs trivial');

subplot(2,3,2);
UQMaT.Plotter.DOS.spin(w,D_up,D_dn,'ax',gca,'title','Spin-resolved DOS (topo)');

subplot(2,3,3);
UQMaT.Plotter.Bands.plot1D(kn,Ek,'ax',gca,'title','Band structure (twisted BC)');

subplot(2,3,4);
UQMaT.Plotter.LDOS.map1D(w,L_t,'ax',gca,'title','LDOS — topological phase');

subplot(2,3,5);
if ~isempty(rho)
    plot(x, rho','LineWidth',1.5);
    xlabel('Site'); ylabel('|\psi|^2');
    title(sprintf('MZM density (E=%.4f)', ev_t.E(idx_t(1))));
    grid on;
end

subplot(2,3,6);
[mp_prof, xp] = mp_t.profile();
plot(xp, mp_prof,'b','LineWidth',1.5);
xlabel('Site'); ylabel('|MP(i)|');
title('Majorana polarization profile'); grid on;
