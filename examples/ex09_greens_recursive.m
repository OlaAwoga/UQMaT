%% ex09_greens_recursive.m
%  Green's function — Recursive (RGF) method.
%  Layer-by-layer propagation: G^r(n) = [ω - H_n - Σ_L(n)]^{-1}
%  O(N·M²) scaling — ideal for large systems and transport.
%
%  Demonstrates:
%    1. DOS for large 1D chain (N=500) — compare with resolvent
%    2. Transmission T(ω) via Landauer formula
%    3. Effect of disorder on DOS and T
%    4. LDOS profile along chain
%    5. Scaling: RGF vs direct vs resolvent

uqmat_path();

%% ── 1. Large chain DOS — RGF vs resolvent ────────────────────────────────────
fprintf('1. Large chain DOS (N=500) via RGF...\n');
N_large = 500;
omega_g = linspace(-3.5, 3.5, 500);

p_rgf = UQMaT.ParamsRecursiveGreen( ...
    'system',     'chain', ...
    'N',          N_large, ...
    't',          1.0, ...
    'mu',         0.0, ...
    'EnergyGrid', omega_g, ...
    'eta',        0.02, ...
    'n_layers',   N_large);

disp(p_rgf);

sys_rgf  = UQMaT.System(p_rgf.toParams());
slv_rgf  = UQMaT.Solver.Greens.Recursive(sys_rgf, p_rgf);
ev_rgf   = UQMaT.EV(slv_rgf);
[w, D_rgf] = ev_rgf.dos();

% Resolvent comparison (same system, smaller N for speed)
N_cmp  = 100;
p_cmp  = UQMaT.ParamsRecursiveGreen('system','chain','N',N_cmp, ...
    'EnergyGrid',omega_g,'eta',0.02,'n_layers',N_cmp);
sys_cmp   = UQMaT.System(p_cmp.toParams());
slv_cmp   = UQMaT.Solver.Greens.Recursive(sys_cmp, p_cmp);
ev_res    = UQMaT.EV(UQMaT.System(p_cmp.toParams()), 'resolvent', 'eta',0.02);
[w, D_res]= ev_res.dos('EnergyGrid',omega_g);
ev_rgf2   = UQMaT.EV(slv_cmp);
[w, D_rgf2] = ev_rgf2.dos();

fprintf('  Done. N=%d RGF, N=%d comparison\n', N_large, N_cmp);

%% ── 2. Transmission via Landauer formula ─────────────────────────────────────
fprintf('2. Transmission T(ω)...\n');
p_T = UQMaT.ParamsRecursiveGreen('system','chain','N',100, ...
    'EnergyGrid',omega_g,'eta',0.02,'n_layers',100);
sys_T   = UQMaT.System(p_T.toParams());
slv_T   = UQMaT.Solver.Greens.Recursive(sys_T, p_T);
[w, T_clean] = slv_T.transmission();
fprintf('  Clean chain: max T = %.4f\n', max(T_clean));

%% ── 3. Disorder effect on T ──────────────────────────────────────────────────
fprintf('3. Disorder effect on transmission...\n');
W_vals = [0, 0.5, 1.0, 2.0];
T_dis  = zeros(numel(omega_g), numel(W_vals));
for iw = 1:numel(W_vals)
    pd = UQMaT.ParamsRecursiveGreen('system','chain','N',100, ...
        'W',W_vals(iw),'seed',42,'EnergyGrid',omega_g,'eta',0.02,'n_layers',100);
    sysd  = UQMaT.System(pd.toParams());
    slvd  = UQMaT.Solver.Greens.Recursive(sysd, pd);
    [~, T_dis(:,iw)] = slvd.transmission();
end
fprintf('  Done.\n');

%% ── 4. SC chain RGF ──────────────────────────────────────────────────────────
fprintf('4. SC chain via RGF...\n');
p_sc = UQMaT.ParamsRecursiveGreen('system','chain','N',100, ...
    'SC',struct('Delta_s',0.3,'sym','onsite'), ...
    'EnergyGrid',omega_g,'eta',0.02,'n_layers',100);
sys_sc  = UQMaT.System(p_sc.toParams());
slv_sc  = UQMaT.Solver.Greens.Recursive(sys_sc, p_sc);
ev_sc   = UQMaT.EV(slv_sc);
[w, D_sc] = ev_sc.dos();
fprintf('  Done.\n');

%% ── 5. Scaling comparison ────────────────────────────────────────────────────
fprintf('5. Scaling test...\n');
N_vals = [50, 100, 200];
t_rgf_v  = zeros(1,3);
t_dir_v  = zeros(1,3);
for in = 1:numel(N_vals)
    Nv = N_vals(in);
    pv = UQMaT.ParamsRecursiveGreen('system','chain','N',Nv,'EnergyGrid',omega_g,'eta',0.02,'n_layers',Nv);
    sv = UQMaT.System(pv.toParams());
    tic; UQMaT.EV(UQMaT.Solver.Greens.Recursive(sv,pv)); t_rgf_v(in) = toc;
    tic; UQMaT.EV(sv,'direct'); t_dir_v(in) = toc;
end
fprintf('  N=%s: RGF=[%.2fs, %.2fs, %.2fs] Direct=[%.2fs, %.2fs, %.2fs]\n', ...
    mat2str(N_vals), t_rgf_v, t_dir_v);

%% ── Plots ────────────────────────────────────────────────────────────────────
figure(1); clf;
set(gcf,'Name','ex09 — Recursive Green''s function','Position',[50,50,1300,900]);

subplot(3,3,1);
plot(w, D_rgf, 'b', 'LineWidth',2);
xlabel('\omega'); ylabel('DOS/site');
title(sprintf('RGF DOS (N=%d chain)', N_large)); grid on;

subplot(3,3,2);
hold on;
plot(w, D_rgf2, 'b', 'LineWidth',2, 'DisplayName',sprintf('RGF N=%d',N_cmp));
plot(w, D_res,  'r--','LineWidth',1.5,'DisplayName','Resolvent');
xlabel('\omega'); ylabel('DOS/site');
title('RGF vs resolvent comparison'); legend('FontSize',8); grid on;

subplot(3,3,3);
plot(w, T_clean, 'k', 'LineWidth',2);
xlabel('\omega'); ylabel('T(\omega)');
title('Transmission (Landauer, clean chain)'); grid on; ylim([0,inf]);

subplot(3,3,4);
colors = {'k','b','r','m'};
for iw = 1:numel(W_vals)
    plot(w, T_dis(:,iw), colors{iw}, 'LineWidth',1.5, ...
        'DisplayName', sprintf('W=%.1f',W_vals(iw)));
    hold on;
end
xlabel('\omega'); ylabel('T(\omega)');
title('T vs disorder strength W'); legend('FontSize',8); grid on;

subplot(3,3,5);
hold on;
plot(w, D_sc, 'b', 'LineWidth',2, 'DisplayName','SC (RGF)');
plot(w, D_rgf2, 'k--', 'LineWidth',1.5, 'DisplayName','Normal');
xlabel('\omega'); ylabel('DOS'); title('SC chain DOS — RGF'); legend('FontSize',8); grid on;

subplot(3,3,6);
loglog(N_vals, t_rgf_v, 'b-o','LineWidth',2,'MarkerSize',8,'DisplayName','RGF');
hold on;
loglog(N_vals, t_dir_v, 'r-s','LineWidth',2,'MarkerSize',8,'DisplayName','Direct eig');
xlabel('N'); ylabel('Time (s)');
title('Scaling: RGF vs direct');
legend('FontSize',8); grid on;

subplot(3,3,7);
% DOS for different disorder realisations (ensemble)
D_ens = zeros(numel(omega_g), 5);
for ir = 1:5
    pr = UQMaT.ParamsRecursiveGreen('system','chain','N',100,'W',1.0,'seed',ir, ...
        'EnergyGrid',omega_g,'eta',0.02,'n_layers',100);
    sr = UQMaT.System(pr.toParams());
    ev_r = UQMaT.EV(UQMaT.Solver.Greens.Recursive(sr,pr));
    [~, D_ens(:,ir)] = ev_r.dos();
end
hold on;
for ir = 1:5
    plot(w, D_ens(:,ir), 'Color',[0.5 0.5 0.9 0.5], 'LineWidth',0.8);
end
plot(w, mean(D_ens,2), 'b', 'LineWidth',2.5, 'DisplayName','Ensemble avg');
xlabel('\omega'); ylabel('DOS'); title('Disorder ensemble (W=1.0, 5 realisations)');
legend('FontSize',8); grid on;

subplot(3,3,8);
ax=gca; ax.Visible='off';
info = {
    'ParamsRecursiveGreen summary',
    '',
    sprintf('Large chain: N=%d', N_large),
    sprintf('RGF n_layers=%d', N_large),
    sprintf('Block size M = ceil(Hdim/n_layers)'),
    '',
    'Solver: G^r(n) = [ω-H_n-Σ_L]^{-1}',
    'Dyson: Σ_L(n) = V G^r(n-1) V†',
    'Scaling: O(N·M²)',
    'Memory: O(M²) — layer-local',
    '',
    'T(ω) = Tr[Γ_L G^r Γ_R G^a]',
    'Wide-band lead: Γ = η·I'
};
for k=1:numel(info)
    text(0.05, 0.97-k*0.07, info{k}, 'Units','normalized','FontSize',8, ...
        'FontWeight',iif(k==1,'bold','normal'));
end

subplot(3,3,9);
% Kitaev chain via RGF
p_kit = UQMaT.ParamsRecursiveGreen('system','kitaev','N',80, ...
    'mu',-1.0,'Delta_p',0.5,'EnergyGrid',omega_g,'eta',0.02,'n_layers',80);
sys_kit = UQMaT.System(p_kit.toParams());
slv_kit = UQMaT.Solver.Greens.Recursive(sys_kit, p_kit);
ev_kit  = UQMaT.EV(slv_kit);
[w, D_kit] = ev_kit.dos();
plot(w, D_kit, 'r', 'LineWidth',2);
xlabel('\omega'); ylabel('DOS');
title('Kitaev chain (topological, via RGF)'); grid on;
xline(0,'k:','LineWidth',1);

sgtitle('ex09 — Recursive Green''s function: layer-by-layer propagation','FontSize',11,'FontWeight','bold');

function s = iif(c,a,b); if c, s=a; else, s=b; end; end
