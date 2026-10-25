%% ex10_greens_kkr.m
%  Green's function — KKR (Korringa-Kohn-Rostoker) method.
%  G^{KKR} = [t(ω)^{-1} - B]^{-1}  via t-matrix and structure constants.
%
%  Demonstrates:
%    1. DOS via KKR for 1D chain and 2D square lattice
%    2. LDOS at impurity site — Friedel oscillations
%    3. CPA self-energy for disordered system
%    4. Compare KKR vs resolvent
%    5. SC gap from KKR Green's function

uqmat_path();

%% ── 1. Chain DOS via KKR ─────────────────────────────────────────────────────
fprintf('1. Chain DOS via KKR...\n');
omega_g = linspace(-3.5, 3.5, 500);

p_kkr = UQMaT.ParamsKKRGreen( ...
    'system',     'chain', ...
    'N',          80, ...
    't',          1.0, ...
    'mu',         0.0, ...
    'EnergyGrid', omega_g, ...
    'eta',        0.02);

disp(p_kkr);

sys_kkr  = UQMaT.System(p_kkr.toParams());
slv_kkr  = UQMaT.Solver.Greens.KKR(sys_kkr, p_kkr);
ev_kkr   = UQMaT.EV(slv_kkr);
[w, D_kkr] = ev_kkr.dos();
fprintf('  Done.\n');

%% ── 2. 2D square via KKR ─────────────────────────────────────────────────────
fprintf('2. 2D square lattice DOS via KKR...\n');
p_sq = UQMaT.ParamsKKRGreen('system','square','Nx',10,'Ny',10, ...
    'EnergyGrid',omega_g,'eta',0.04);
sys_sq   = UQMaT.System(p_sq.toParams());
slv_sq   = UQMaT.Solver.Greens.KKR(sys_sq, p_sq);
ev_sq    = UQMaT.EV(slv_sq);
[w, D_sq] = ev_sq.dos();

% Reference: resolvent
ev_sq_res = UQMaT.EV(sys_sq, 'resolvent', 'eta',0.04);
[w, D_sq_res] = ev_sq_res.dos('EnergyGrid',omega_g);
fprintf('  Done.\n');

%% ── 3. LDOS via KKR at multiple sites ───────────────────────────────────────
fprintf('3. Site-resolved LDOS via KKR...\n');
sites_ldos = [1, 20, 40, 80];
[w, L_kkr] = slv_kkr.ldos(sites_ldos);
fprintf('  Done.\n');

%% ── 4. Disorder DOS — Anderson vs CPA ───────────────────────────────────────
fprintf('4. CPA self-energy for disordered system...\n');
p_cpa = UQMaT.ParamsKKRGreen('system','chain','N',80, ...
    'W',1.0,'seed',42,'EnergyGrid',omega_g,'eta',0.02, ...
    'use_cpa',true,'cpa_tol',1e-6,'cpa_maxit',100);
sys_cpa = UQMaT.System(p_cpa.toParams());
slv_cpa = UQMaT.Solver.Greens.KKR(sys_cpa, p_cpa);
ev_cpa  = UQMaT.EV(slv_cpa);
[w, D_cpa]   = ev_cpa.dos();
[w, Sig_cpa] = slv_cpa.self_energy();

% Disordered system via direct diag (ensemble average)
D_dis_ens = zeros(numel(omega_g), 5);
for ir = 1:5
    pd  = UQMaT.ParamsKKRGreen('system','chain','N',80,'W',1.0,'seed',ir, ...
        'EnergyGrid',omega_g,'eta',0.02);
    evd = UQMaT.EV(UQMaT.System(pd.toParams()),'direct');
    [~, D_dis_ens(:,ir)] = evd.dos('EnergyGrid',omega_g,'eta',0.02);
end
D_dis_avg = mean(D_dis_ens, 2);
fprintf('  CPA done.\n');

%% ── 5. SC chain via KKR ──────────────────────────────────────────────────────
fprintf('5. SC chain via KKR...\n');
p_sc = UQMaT.ParamsKKRGreen('system','chain','N',80, ...
    'SC',struct('Delta_s',0.3,'sym','onsite'),'EnergyGrid',omega_g,'eta',0.02);
sys_sc  = UQMaT.System(p_sc.toParams());
slv_sc  = UQMaT.Solver.Greens.KKR(sys_sc, p_sc);
ev_sc   = UQMaT.EV(slv_sc);
[w, D_sc] = ev_sc.dos();
fprintf('  Done.\n');

%% ── Plots ────────────────────────────────────────────────────────────────────
figure(1); clf;
set(gcf,'Name','ex10 — KKR Green''s function','Position',[50,50,1300,900]);

subplot(3,3,1);
plot(w, D_kkr, 'b', 'LineWidth',2);
xlabel('\omega'); ylabel('DOS/site');
title('1D chain DOS — KKR method'); grid on;

subplot(3,3,2);
hold on;
plot(w, D_sq,     'b',  'LineWidth',2,   'DisplayName','KKR');
plot(w, D_sq_res, 'r--','LineWidth',1.5, 'DisplayName','Resolvent');
xlabel('\omega'); ylabel('DOS');
title('2D square: KKR vs resolvent'); legend('FontSize',8); grid on;

subplot(3,3,3);
for si = 1:size(L_kkr,2)
    plot(w, L_kkr(:,si), 'LineWidth',1.5, ...
        'DisplayName', sprintf('site %d',sites_ldos(si)));
    hold on;
end
xlabel('\omega'); ylabel('LDOS');
title('Site-resolved LDOS — KKR'); legend('FontSize',7); grid on;

subplot(3,3,4);
hold on;
plot(w, D_cpa,    'b',  'LineWidth',2,   'DisplayName','KKR-CPA');
plot(w, D_dis_avg,'r--','LineWidth',1.5, 'DisplayName','Ensemble avg (direct)');
xlabel('\omega'); ylabel('DOS');
title('Disordered chain: CPA vs ensemble'); legend('FontSize',8); grid on;

subplot(3,3,5);
hold on;
plot(w, real(Sig_cpa), 'b', 'LineWidth',1.5, 'DisplayName','Re \Sigma');
plot(w, imag(Sig_cpa), 'r', 'LineWidth',1.5, 'DisplayName','Im \Sigma');
xlabel('\omega'); ylabel('\Sigma(\omega)');
title('CPA self-energy \Sigma(ω)'); legend('FontSize',8); grid on;

subplot(3,3,6);
hold on;
plot(w, D_sc,  'b', 'LineWidth',2, 'DisplayName','SC');
plot(w, D_kkr, 'k--','LineWidth',1.5,'DisplayName','Normal');
xlabel('\omega'); ylabel('DOS');
title('SC gap in KKR Green''s function'); legend('FontSize',8); grid on;

subplot(3,3,7);
% Impurity effect: single magnetic impurity at site 40
p_imp = UQMaT.ParamsKKRGreen('system','chain','N',80,'imp_pos',40, ...
    'imp_potential',2.0,'EnergyGrid',omega_g,'eta',0.02);
sys_imp = UQMaT.System(p_imp.toParams());
slv_imp = UQMaT.Solver.Greens.KKR(sys_imp, p_imp);
[w, L_imp] = slv_imp.ldos([38,39,40,41,42]);
imagesc([38,39,40,41,42], w, L_imp);
axis xy; colorbar; colormap(gca,hot);
xlabel('Site'); ylabel('\omega');
title('LDOS near impurity (site 40)');
caxis([0, quantile(L_imp(:),0.99)]);

subplot(3,3,8);
% Friedel oscillations: |G(40,j)| vs j
H_imp  = full(sys_imp.H);
z_0    = 0.2 + 1i*0.02;
G_imp  = (z_0*eye(size(H_imp)) - H_imp) \ eye(size(H_imp));
plot(1:80, abs(G_imp(40,:)), 'b.-', 'LineWidth',1.5);
xlabel('Site j'); ylabel('|G(40,j,\omega_0)|');
title('Friedel oscillations around impurity');
set(gca,'YScale','log'); grid on;

subplot(3,3,9);
ax=gca; ax.Visible='off';
info = {
    'ParamsKKRGreen summary',
    '',
    'G^{KKR} = [t(ω)^{-1} - B]^{-1}',
    't = single-site t-matrix',
    'B = off-diagonal (structure const.)',
    '',
    'CPA self-energy:',
    'Σ s.t. <G_{CPA}> = G_0(ω-Σ)',
    'for configurationally averaged G',
    '',
    sprintf('N=%d chain  eta=%.3f', p_kkr.N, p_kkr.eta),
    sprintf('Grid: %d pts  [%.2f,%.2f]', numel(omega_g), omega_g(1), omega_g(end))
};
for k=1:numel(info)
    text(0.05, 0.97-k*0.075, info{k}, 'Units','normalized','FontSize',8, ...
        'FontWeight',iif(k==1,'bold','normal'));
end

sgtitle('ex10 — KKR Green''s function: t-matrix + structure constants','FontSize',11,'FontWeight','bold');

function s = iif(c,a,b); if c, s=a; else, s=b; end; end
