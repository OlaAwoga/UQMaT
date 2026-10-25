%% ex04_kitaev_chain.m
%  Kitaev p-wave SC: Majorana zero modes, Majorana polarization, phase diagram.
%  Demonstrates: single EV diagonalisation feeding DOS, MP, and topological invariants.

%
%  Solver options (change 'direct' to any of these):
%    ev = UQMaT.EV(sys, 'direct')                        % dense eig()   — default, all states
%    ev = UQMaT.EV(sys, 'lanczos', 'k',20, 'sigma',0)   % sparse eigs() — k states near E=0
%    ev = UQMaT.EV(sys, 'resolvent', 'eta',0.02)         % G(ω) via inversion
%    ev = UQMaT.EV(sys, 'chebyshev', 'M',512)            % KPM — very large systems
uqmat_path();

N = 80;
t = 1.0;

% ---- Topological (|mu|<2t) — single diagonalisation ----
p_topo = UQMaT.Params('system','kitaev','N',N,'t',t,'mu',-1.0,'Delta_p',0.5,'bc','obc');
sys_t = UQMaT.System(p_topo);
ev_t  = UQMaT.EV(sys_t, 'direct');   % dense eig(), all states

% All observables read from the same ev_t — eig() called once
[w, D_topo]  = ev_t.dos('range',[-3,3]);
idx          = ev_t.edge_idx();
fprintf('In-gap states: %d  E = %s\n', numel(idx), num2str(ev_t.E(idx)','%.5f '));

% Majorana polarization
mp = ev_t.mp();
disp(mp)

C                     = mp.charge(idx);
MP                    = mp.local(idx);
[MP_L, MP_R, MP_LR]  = mp.left_right(idx);
MP_nl                 = mp.nonlocal(idx);
[NMP, MP_LR_nl]       = mp.nonlocal_total(idx);
[mp_prof, xp]         = mp.profile();
[nl_prof, xnl]        = mp.nonlocal_profile(N, idx);

fprintf('MP_L=%.4f  MP_R=%.4f  MP_LR=%.4f  (< 0 → topological)\n',...
    MP_L(1), MP_R(1), MP_LR(1));
fprintf('NMP=%.4f  MP_LR_nl=%.4f\n', NMP, MP_LR_nl);

% ---- Trivial (|mu|>2t) ----
p_triv = UQMaT.Params('system','kitaev','N',N,'t',t,'mu',-2.5,'Delta_p',0.5,'bc','obc');
sys_v = UQMaT.System(p_triv);
ev_v  = UQMaT.EV(sys_v, 'direct');   % dense eig(), all states
[w, D_triv] = ev_v.dos('range',[-3,3]);

% ---- Bulk invariants ----
inv_t = UQMaT.Observables.Topological.Bulk.Invariants(...
    UQMaT.System(UQMaT.Params('system','kitaev','N',N,'t',t,'mu',-1.0,...
    'Delta_p',0.5,'bc','pbc')));
disp(inv_t);

% ---- Real-space (from cached ev_t) ----
rs_t     = ev_t.realspace_invariants();
[rho, x] = rs_t.edge_density();
[ipr_v, E_v] = rs_t.ipr();

% ---- Band structure ----
[kn, Ek] = UQMaT.Models.Chains1D.ChainBC.bands(...
    UQMaT.System(UQMaT.Params('system','kitaev','N',1,'t',t,'mu',-1.0,'Delta_p',0.5)), 300);

% ---- Phase diagram ----
p0 = UQMaT.Params('system','kitaev','N',30,'t',t,'bc','pbc');
pd = UQMaT.Observables.Topological.PhaseDiagram(p0,'Npts',30);
[mu_v, Dp_v, nu_map, gap_map] = pd.kitaev_phase([-3 3],[0.05 1.5]);

% ---- Plot ----
figure(1); clf;
set(gcf,'Name','ex04 — Kitaev chain + Majorana polarization');

subplot(2,4,1);
UQMaT.Plotter.DOS.compare(w,{D_topo,D_triv},...
    'labels',{'Topo |μ|<2t','Trivial |μ|>2t'},'fill',true,'ax',gca,...
    'title','DOS: topological vs trivial');

subplot(2,4,2);
UQMaT.Plotter.Bands.plot1D(kn,Ek,'ax',gca,'title','BdG bands');

subplot(2,4,3);
if ~isempty(rho)
    plot(x, rho','LineWidth',2);
    xlabel('Site'); ylabel('|\psi|^2');
    title(sprintf('MZM density  E≈%.5f', ev_t.E(idx(1))));
    grid on;
end

subplot(2,4,4);
plot(xp, mp_prof,'b','LineWidth',2);
xlabel('Site'); ylabel('|MP(i)|');
title('Majorana polarization |MP(i)|'); grid on;

subplot(2,4,5);
sites_v = 1:N;
imagesc(sites_v, sites_v, MP_nl);
axis xy equal tight; colorbar; colormap(gca,redblue_or_parula());
xlabel('j'); ylabel('i');
title('Nonlocal MP_{nl}(i,j)');

subplot(2,4,6);
plot(xnl, nl_prof,'r','LineWidth',2);
xlabel('Site i'); ylabel('MP_{nl}(i,N)');
title(sprintf('Nonlocal profile (j_{ref}=%d)',N)); grid on;

subplot(2,4,7);
UQMaT.Plotter.Phase.winding_map(mu_v, Dp_v, nu_map,'ax',gca,...
    'xlabel','\Delta_p','ylabel','\mu','title','Winding number \nu');

subplot(2,4,8);
UQMaT.Plotter.Phase.gap_map(mu_v, Dp_v, gap_map,'ax',gca,...
    'xlabel','\Delta_p','ylabel','\mu','title','Quasiparticle gap');

function c = redblue_or_parula()
    try, c = redblue(); catch, c = parula(); end
end
