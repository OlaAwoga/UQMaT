%% ex08_greens_resolvent.m
%  Green's function — Resolvent (direct inversion) method.
%  G^r(ω) = (ω + iη - H)^{-1}  exact, O(N³) per ω point.
%
%  Demonstrates:
%    1. DOS via Tr Im G — compare with eigenvalue method
%    2. LDOS at specific sites via G(i,i,ω)
%    3. Effect of broadening η on spectral resolution
%    4. SC gap from G — particle-hole structure
%    5. Off-diagonal G(i,j) — spatial coherence length

uqmat_path();

%% ── Parameters via ParamsResolventGreen ──────────────────────────────────────
fprintf('Building resolvent Green''s function solver...\n');

p_rg = UQMaT.ParamsResolventGreen( ...
    'system',     'chain', ...
    'N',          80, ...
    't',          1.0, ...
    'mu',         0.0, ...
    'EnergyGrid', linspace(-3.5, 3.5, 600), ...
    'eta',        0.02);

disp(p_rg);

sys_rg = UQMaT.System(p_rg.toParams());
ev_rg  = UQMaT.EV(sys_rg, 'resolvent', 'eta', p_rg.eta);

[w, D_rg] = ev_rg.dos();

% Compare with direct diagonalisation
ev_dir = UQMaT.EV(sys_rg, 'direct');
[w, D_dir] = ev_dir.dos('EnergyGrid', p_rg.EnergyGrid, 'eta', p_rg.eta);

fprintf('  Normal chain: resolvent vs direct\n');

%% ── SC chain via resolvent ───────────────────────────────────────────────────
p_sc = UQMaT.ParamsResolventGreen( ...
    'system',     'chain', ...
    'N',          80, ...
    'SC',         struct('Delta_s',0.3,'sym','onsite'), ...
    'EnergyGrid', linspace(-3.5, 3.5, 600), ...
    'eta',        0.02);

sys_sc = UQMaT.System(p_sc.toParams());
ev_sc  = UQMaT.EV(sys_sc, 'resolvent', 'eta', p_sc.eta);

[w, D_sc]     = ev_sc.dos();
[w, D_sc_e, D_sc_h] = ev_sc.dos_ph();
[w, L_sc]     = ev_sc.ldos([1, 20, 40, 80]);  % sites at edge and bulk

fprintf('  SC chain resolved\n');

%% ── Effect of η on spectral resolution ──────────────────────────────────────
eta_vals = [0.005, 0.02, 0.05, 0.1];
D_eta    = zeros(numel(w), numel(eta_vals));
for ie = 1:numel(eta_vals)
    ev_eta = UQMaT.EV(sys_rg, 'resolvent', 'eta', eta_vals(ie));
    [~, D_eta(:,ie)] = ev_eta.dos('EnergyGrid', p_rg.EnergyGrid);
end
fprintf('  η comparison done\n');

%% ── Off-diagonal G(1,j) — coherence length ───────────────────────────────────
% G(1,j,ω) = coherence between site 1 and site j at energy ω
H_full = full(sys_rg.H);
omega0 = 0.5;    % energy near band edge
eta0   = 0.02;
z      = omega0 + 1i*eta0;
G_full = (z*eye(size(H_full)) - H_full) \ eye(size(H_full));
G_offdiag = abs(G_full(1,:));   % |G(1,j)| for all j

fprintf('  Off-diagonal G computed\n');

%% ── Plots ────────────────────────────────────────────────────────────────────
figure(1); clf;
set(gcf,'Name','ex08 — Resolvent Green''s function','Position',[50,50,1300,900]);

subplot(3,3,1);
hold on;
plot(w, D_rg,  'b', 'LineWidth',2, 'DisplayName','Resolvent G');
plot(w, D_dir, 'r--','LineWidth',1.5,'DisplayName','Direct eig');
xlabel('\omega'); ylabel('DOS/site');
title('Resolvent vs direct diagonalisation'); legend('FontSize',8); grid on;

subplot(3,3,2);
hold on;
plot(w, D_sc,   'k',  'LineWidth',2, 'DisplayName','Total');
plot(w, D_sc_e, 'b',  'LineWidth',1.5,'DisplayName','Electron');
plot(w, D_sc_h, 'r--','LineWidth',1.5,'DisplayName','Hole');
xlabel('\omega'); ylabel('DOS'); title('SC chain — particle-hole DOS');
legend('FontSize',8); grid on;

subplot(3,3,3);
for ie = 1:numel(eta_vals)
    plot(w, D_eta(:,ie)/max(D_eta(:,ie)+eps), 'LineWidth',1.5, ...
        'DisplayName', sprintf('\\eta=%.3f',eta_vals(ie)));
    hold on;
end
xlabel('\omega'); ylabel('DOS (norm.)');
title('Effect of broadening \eta'); legend('FontSize',7); grid on;

subplot(3,3,4);
imagesc(1:4, w, L_sc);
axis xy; colorbar; colormap(gca,hot);
xlabel('Site'); ylabel('\omega');
title('LDOS at sites 1, 20, 40, 80 (SC chain)');
set(gca,'XTick',1:4,'XTickLabel',{'1','20','40','80'});
caxis([0, quantile(L_sc(:),0.99)]);

subplot(3,3,5);
plot(1:p_rg.N, G_offdiag(1:p_rg.N), 'b.-', 'LineWidth',1.5);
xlabel('Site j'); ylabel('|G(1,j,\omega_0)|');
title(sprintf('Off-diagonal coherence (\\omega_0=%.2f)', omega0));
set(gca,'YScale','log'); grid on;

subplot(3,3,6);
% DOS on 2D square via resolvent
p2d = UQMaT.ParamsResolventGreen('system','square','Nx',10,'Ny',10, ...
    'EnergyGrid',linspace(-5,5,500),'eta',0.05);
ev2d = UQMaT.EV(UQMaT.System(p2d.toParams()), 'resolvent', 'eta', p2d.eta);
[w2, D2d] = ev2d.dos();
plot(w2, D2d, 'k', 'LineWidth',2);
xlabel('\omega'); ylabel('DOS/site');
title('2D square lattice DOS (resolvent)');
grid on; xline(0,'r:');

subplot(3,3,7);
% Resolvent imaginary part = spectral function
w_grid = linspace(-3,3,200);
G_diag_all = zeros(numel(w_grid), p_rg.N);
for iw = 1:numel(w_grid)
    z_  = w_grid(iw) + 1i*0.02;
    Gz  = (z_*eye(size(H_full)) - H_full) \ eye(size(H_full));
    G_diag_all(iw,:) = -imag(diag(Gz))/pi;
end
imagesc(1:p_rg.N, w_grid, G_diag_all);
axis xy; colorbar; colormap(gca,hot);
xlabel('Site i'); ylabel('\omega');
title('A(i,i,\omega) = -ImG_{ii}/\pi (spectral function)');
caxis([0, quantile(G_diag_all(:),0.99)]);

subplot(3,3,8);
% Real vs imaginary part of G(ω)
plot(w, real(D_rg*pi), 'b', 'LineWidth',1.5, 'DisplayName','Re G (× \pi)');
hold on;
plot(w, D_rg, 'r', 'LineWidth',1.5, 'DisplayName','-Im G/\pi = DOS');
xlabel('\omega'); ylabel('G(\omega)');
title('Real and imaginary parts of Tr G');
legend('FontSize',8); grid on;

subplot(3,3,9);
ax = gca;
text(0.05,0.90,sprintf('ParamsResolventGreen summary'),'FontSize',9,'FontWeight','bold','Units','normalized');
info = {
    sprintf('System: chain, N=%d', p_rg.N),
    sprintf('t=%.2f  mu=%.2f', p_rg.t, p_rg.mu),
    sprintf('Energy grid: %d points', numel(p_rg.EnergyGrid)),
    sprintf('[%.2f, %.2f]', p_rg.EnergyGrid(1), p_rg.EnergyGrid(end)),
    sprintf('eta=%.4f', p_rg.eta),
    '',
    'Solver: Resolvent',
    'G(ω) = (ω+iη-H)^{-1}',
    'Exact — O(N³) per ω',
    'Best for N ≲ 1000'
};
for k=1:numel(info)
    text(0.05, 0.82-k*0.08, info{k}, 'Units','normalized','FontSize',8);
end
ax.Visible = 'off';

sgtitle('ex08 — Resolvent Green''s function: G(ω) = (ω+iη−H)^{−1}','FontSize',11,'FontWeight','bold');
