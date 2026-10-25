%% ex03_ssh_chain.m
%  Su-Schrieffer-Heeger chain: trivial vs topological phases.
%  Demonstrates: Zak phase, edge states, band structure, phase diagram.

%
%  Solver options (change 'direct' to any of these):
%    ev = UQMaT.EV(sys, 'direct')                        % dense eig()   — default, all states
%    ev = UQMaT.EV(sys, 'lanczos', 'k',20, 'sigma',0)   % sparse eigs() — k states near E=0
%    ev = UQMaT.EV(sys, 'resolvent', 'eta',0.02)         % G(ω) via inversion
%    ev = UQMaT.EV(sys, 'chebyshev', 'M',512)            % KPM — very large systems
uqmat_path();

N = 40;

% ---- Trivial and topological — single diagonalisation each ----
p_triv = UQMaT.Params('system','ssh','N',N,'v',1.0,'w',0.3,'bc','obc');
p_topo = UQMaT.Params('system','ssh','N',N,'v',0.3,'w',1.0,'bc','obc');

sys_triv = UQMaT.System(p_triv);
ev_triv  = UQMaT.EV(sys_triv, 'direct');  % dense eig(), all states
sys_topo = UQMaT.System(p_topo);
ev_topo  = UQMaT.EV(sys_topo, 'direct');  % dense eig(), all states

[w, D_triv] = ev_triv.dos('range',[-2,2]);
[w, D_topo] = ev_topo.dos('range',[-2,2]);

% ---- Bulk invariants (PBC) ----
inv_triv = UQMaT.Observables.Topological.Bulk.Invariants(...
    UQMaT.System(UQMaT.Params('system','ssh','N',N,'v',1.0,'w',0.3,'bc','pbc')));
inv_topo = UQMaT.Observables.Topological.Bulk.Invariants(...
    UQMaT.System(UQMaT.Params('system','ssh','N',N,'v',0.3,'w',1.0,'bc','pbc')));

phi_triv = inv_triv.zak();
phi_topo = inv_topo.zak();
fprintf('Zak — trivial: %.4f rad   topological: %.4f rad\n', phi_triv, phi_topo);

% ---- Real-space (from cached EV) ----
rs_topo  = ev_topo.realspace_invariants();
[rho, x] = rs_topo.edge_density();

% ---- Band structure ----
[kn, Ek_triv] = UQMaT.Models.Chains1D.ChainBC.bands(...
    UQMaT.System(UQMaT.Params('system','ssh','v',1.0,'w',0.3)), 300);
[kn, Ek_topo] = UQMaT.Models.Chains1D.ChainBC.bands(...
    UQMaT.System(UQMaT.Params('system','ssh','v',0.3,'w',1.0)), 300);

% ---- Phase diagram ----
p0 = UQMaT.Params('system','ssh','N',20,'bc','pbc');
pd = UQMaT.Observables.Topological.PhaseDiagram(p0,'Npts',30);
[v_v, w_v, phi_map, gap_map] = pd.ssh_phase([0.1 2],[0.1 2]);

% ---- Plot ----
figure(1); clf;
set(gcf,'Name','ex03 — SSH chain');

subplot(2,3,1);
UQMaT.Plotter.Spectrum.compare({UQMaT.System(p_triv), UQMaT.System(p_topo)},...
    'labels',{'Trivial (v>w)','Topological (w>v)'},'ax',gca,'title','Spectra');

subplot(2,3,2);
UQMaT.Plotter.DOS.compare(w,{D_triv,D_topo},...
    'labels',{'Trivial','Topological'},'fill',true,'ax',gca,'title','DOS');

subplot(2,3,3);
hold on;
UQMaT.Plotter.Bands.plot1D(kn,Ek_triv,'ax',gca,'color',[0 0.4 0.8],'title','Band structure');
UQMaT.Plotter.Bands.plot1D(kn,Ek_topo,'ax',gca,'color',[0.8 0.2 0]);
legend('Trivial','Topological');

subplot(2,3,4);
if ~isempty(rho)
    plot(x, rho','LineWidth',1.5);
    xlabel('Site'); ylabel('|\psi|^2');
    title('Edge states (topological)'); grid on;
end

subplot(2,3,5);
imagesc(w_v, v_v, mod(phi_map,pi)>pi/2);
axis xy; colormap(gca,[0.9 0.9 0.9; 0.2 0.5 0.9]);
xlabel('w (inter)'); ylabel('v (intra)');
title('Zak phase (blue=\pi, grey=0)'); axis square;

subplot(2,3,6);
UQMaT.Plotter.Phase.gap_map(v_v,w_v,gap_map,'ax',gca,...
    'xlabel','w','ylabel','v','title','Gap map');
hold on; plot(v_v,v_v,'w--','LineWidth',1.5);
