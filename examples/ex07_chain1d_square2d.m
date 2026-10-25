%% ex07_chain1d_square2d.m
%  1D chain coupled to a 2D square lattice substrate.
%
%  Physical scenarios:
%    1. Normal chain on normal substrate  — hybridisation and LDOS
%    2. Normal chain on SC substrate      — proximity-induced gap in wire
%    3. Topological wire (Rashba+Zeeman)  — on SC substrate (MZM in wire)
%    4. Coupling strength Gamma scan      — gap and in-gap states vs Gamma
%    5. Chain position: edge vs bulk      — substrate surface effect
%    6. Majorana polarization             — confirm topological MZM in wire
%
%  Geometry:
%    Substrate : Nx x Ny square lattice (SC, s-wave)
%    Chain     : N_chain sites running along x, at substrate row chain_row
%    Coupling  : Gamma (on-site, chain site k <-> substrate site (k, chain_row))
%
%  Solver options (change 'direct' to any of these):
%    ev = UQMaT.EV(sys, 'direct')                        % dense eig()   — default
%    ev = UQMaT.EV(sys, 'lanczos', 'k',20, 'sigma',0)   % sparse eigs() — k states near E=0
%    ev = UQMaT.EV(sys, 'resolvent', 'eta',0.02)         % G(ω) via inversion
%    ev = UQMaT.EV(sys, 'chebyshev', 'M',512)            % KPM — very large systems

uqmat_path();

%% ── Parameters ───────────────────────────────────────────────────────────────
Nx          = 20;
Ny          = 20;
N_chain     = 20;
t_sub       = 1.0;
t_wire      = 0.8;
mu_sub      = 0.0;
mu_wire     = -0.5;
Delta_s     = 0.35;    % substrate s-wave gap
Gamma       = 0.3;     % chain-substrate coupling
lambda_R    = 0.5;     % Rashba SOC on wire
Bz          = 0.5;     % Zeeman on wire  (topological when Bz > Delta_s)
chain_row_bulk = 10;   % chain in substrate bulk
chain_row_edge = 1;    % chain at substrate edge
eta         = 0.03;
Ns2d        = Nx * Ny;   % number of substrate physical sites

%% ── Scenario 1: Normal chain on normal substrate ────────────────────────────
fprintf('Scenario 1: normal chain on normal substrate...\n');
p1   = ex07_make_params(false, false, chain_row_bulk, Gamma, ...
                        Nx, Ny, N_chain, t_sub, t_wire, mu_sub, mu_wire, ...
                        Delta_s, lambda_R, Bz, eta);
sys1 = UQMaT.System(p1);
ev1  = UQMaT.EV(sys1, 'direct');   % dense eig(), all states

% Wire LDOS row indices depend on Hdim (Ns2d rows for normal substrate,
% 4*Ns2d for BdG substrate). Always use the LAST N_chain rows.
wire_rows1 = (ev1.Hdim - N_chain + 1) : ev1.Hdim;

[w, D1]       = ev1.dos('range',[-3,3]);
[w, L1_chain] = ev1.ldos(wire_rows1, 'range',[-3,3]);
fprintf('  Hdim=%d  wire rows %d:%d\n', ev1.Hdim, wire_rows1(1), wire_rows1(end));

%% ── Scenario 2: Normal chain on SC substrate (proximity effect) ──────────────
fprintf('Scenario 2: chain on SC substrate (proximity effect)...\n');
p2   = ex07_make_params(true, false, chain_row_bulk, Gamma, ...
                        Nx, Ny, N_chain, t_sub, t_wire, mu_sub, mu_wire, ...
                        Delta_s, lambda_R, Bz, eta);
sys2 = UQMaT.System(p2);
ev2  = UQMaT.EV(sys2, 'direct');   % dense eig(), all states

wire_rows2 = (ev2.Hdim - N_chain + 1) : ev2.Hdim;

[w, D2]       = ev2.dos('range',[-3,3]);
[w, L2_chain] = ev2.ldos(wire_rows2, 'range',[-3,3]);
% Note: dos_ph is not called here because the heterostructure has a mixed
% BdG basis (SC substrate + normal chain) — use full dos() instead.
% The pair amplitude is extracted from the substrate portion.
F2            = ev2.pair_amplitude();

% pair_amplitude returns (Ns2d+N_chain)×1 in physical site order
F2_sub    = F2(1:Ns2d);
F2_sub_2D = reshape(abs(F2_sub), Ny, Nx);
if numel(F2) >= Ns2d + N_chain
    F2_wire = abs(F2(Ns2d+1:Ns2d+N_chain));
else
    F2_wire = zeros(N_chain,1);
end

[~,i0] = min(abs(w));
fprintf('  Hdim=%d  wire LDOS at ω=0: %.4f\n', ev2.Hdim, mean(L2_chain(i0,:)));

%% ── Scenario 3: Topological wire (Rashba+Zeeman) on SC substrate ─────────────
fprintf('Scenario 3: topological wire (Rashba+Zeeman) on SC substrate...\n');
p3   = ex07_make_params(true, true, chain_row_bulk, Gamma, ...
                        Nx, Ny, N_chain, t_sub, t_wire, mu_sub, mu_wire, ...
                        Delta_s, lambda_R, Bz, eta);
sys3 = UQMaT.System(p3);
ev3  = UQMaT.EV(sys3, 'direct');   % dense eig(), all states

wire_rows3 = (ev3.Hdim - N_chain + 1) : ev3.Hdim;

[w, D3]       = ev3.dos('range',[-3,3]);
[w, L3_chain] = ev3.ldos(wire_rows3, 'range',[-3,3]);

idx3 = ev3.edge_idx();
fprintf('  Hdim=%d  in-gap states: %d\n', ev3.Hdim, numel(idx3));
if ~isempty(idx3)
    fprintf('  E = '); fprintf('%.5f  ', ev3.E(idx3)); fprintf('\n');
end

%% ── Scenario 4: Coupling strength Gamma scan ─────────────────────────────────
fprintf('Scenario 4: Gamma scan (12 points)...\n');
Gamma_vals       = linspace(0, 1.0, 12);
gap_vs_gamma     = zeros(size(Gamma_vals));
n_ingap_vs_gamma = zeros(size(Gamma_vals));

for ig = 1:numel(Gamma_vals)
    pg   = ex07_make_params(true, true, chain_row_bulk, Gamma_vals(ig), ...
                            Nx, Ny, N_chain, t_sub, t_wire, mu_sub, mu_wire, ...
                            Delta_s, lambda_R, Bz, eta);
    sysg = UQMaT.System(pg);
    evg  = UQMaT.EV(sysg, 'direct');  % dense eig(), all states
    idx_g = evg.edge_idx();
    n_ingap_vs_gamma(ig) = numel(idx_g);
    pos = evg.E(evg.E > 1e-6);
    neg = evg.E(evg.E < -1e-6);
    if ~isempty(pos) && ~isempty(neg)
        gap_vs_gamma(ig) = pos(1) - neg(end);
    end
end
fprintf('  Done.\n');

%% ── Scenario 5: Chain position — edge vs bulk ────────────────────────────────
fprintf('Scenario 5: edge (iy=1) vs bulk (iy=%d) coupling...\n', chain_row_bulk);
p5e  = ex07_make_params(true, true, chain_row_edge, Gamma, ...
                        Nx, Ny, N_chain, t_sub, t_wire, mu_sub, mu_wire, ...
                        Delta_s, lambda_R, Bz, eta);
p5b  = ex07_make_params(true, true, chain_row_bulk, Gamma, ...
                        Nx, Ny, N_chain, t_sub, t_wire, mu_sub, mu_wire, ...
                        Delta_s, lambda_R, Bz, eta);
sys5e = UQMaT.System(p5e);
sys5b = UQMaT.System(p5b);
ev5e  = UQMaT.EV(sys5e, 'direct');   % dense eig(), all states
ev5b  = UQMaT.EV(sys5b, 'direct');   % dense eig(), all states

[w, D5e] = ev5e.dos('range',[-3,3]);
[w, D5b] = ev5b.dos('range',[-3,3]);
fprintf('  Edge in-gap: %d   Bulk in-gap: %d\n', ...
    numel(ev5e.edge_idx()), numel(ev5b.edge_idx()));

%% ── Majorana polarization (scenario 3, wire sites only) ──────────────────────
% MajoranaPolarization is built with N = sys3.Ndim = Ns2d + N_chain.
% left_right() splits at N/2 which bisects the SUBSTRATE, not the wire.
% We therefore compute MP only on wire eigenvector components.
mp3 = ev3.mp();    % MajoranaPolarization for all Ns2d+N_chain physical sites

if ~isempty(idx3)
    % Extract MP profile and take the wire portion (last N_chain sites)
    [mp_prof, ~] = mp3.profile();
    mp_wire = mp_prof(end-N_chain+1:end);   % wire sites only

    % Left/right split ON THE WIRE only: first half vs second half of chain
    MP_wire_mat = mp3.local(idx3);               % Nidx × (Ns2d+N_chain)
    MP_wire_only = MP_wire_mat(:, end-N_chain+1:end);  % Nidx × N_chain
    MP_L_wire = sum(MP_wire_only(:, 1:floor(N_chain/2)),   2);
    MP_R_wire = sum(MP_wire_only(:, floor(N_chain/2)+1:end), 2);
    MP_LR     = MP_L_wire .* MP_R_wire;

    fprintf('Wire MP_L=%.4f  MP_R=%.4f  MP_LR=%.4f  (< 0 = topological MZM)\n', ...
        MP_L_wire(1), MP_R_wire(1), MP_LR(1));
else
    mp_wire = zeros(N_chain,1);
    MP_LR   = 0;
    fprintf('No in-gap states found — trivial or gap too large.\n');
end

%% ── Plots ────────────────────────────────────────────────────────────────────
figure(1); clf;
set(gcf,'Name','ex07 — 1D chain on 2D square lattice','Position',[50,50,1400,900]);

% Row 1 — DOS and LDOS comparisons
subplot(3,4,1);
hold on;
plot(w, D1/max(D1+eps), 'k', 'LineWidth',1.5, 'DisplayName','Normal/Normal');
plot(w, D2/max(D2+eps), 'b', 'LineWidth',1.5, 'DisplayName','Wire/SC sub');
plot(w, D3/max(D3+eps), 'r', 'LineWidth',1.5, 'DisplayName','Topo/SC sub');
xlabel('\omega'); ylabel('DOS (norm.)');
title('DOS: three scenarios'); legend('Location','northwest','FontSize',7); grid on;

subplot(3,4,2);
hold on;
plot(w, mean(L1_chain,2), 'k', 'LineWidth',1.5, 'DisplayName','Normal');
plot(w, mean(L2_chain,2), 'b', 'LineWidth',1.5, 'DisplayName','Prox. SC');
plot(w, mean(L3_chain,2), 'r', 'LineWidth',1.5, 'DisplayName','Topo (Rashba+B)');
xline(0,'k:','LineWidth',0.8);
xlabel('\omega'); ylabel('LDOS/site');
title('Wire LDOS (avg over chain)'); legend('FontSize',7); grid on;

subplot(3,4,3);
hold on;
plot(w, D2/max(D2+eps), 'b',  'LineWidth',1.5, 'DisplayName','Prox. SC');
plot(w, D1/max(D1+eps), 'k--','LineWidth',1.2, 'DisplayName','Normal');
xlabel('\omega'); ylabel('DOS (norm.)');
title('SC substrate proximity effect on DOS'); legend('FontSize',8); grid on;
% Note: particle-hole resolved LDOS requires a uniform BdG system.
% For the heterostructure use ev.dos() on the full system.

subplot(3,4,4);
hold on;
plot(w, D5e/max(D5e+eps), 'b', 'LineWidth',1.5, 'DisplayName','Edge (iy=1)');
plot(w, D5b/max(D5b+eps), 'r', 'LineWidth',1.5, 'DisplayName','Bulk (iy=10)');
xlabel('\omega'); ylabel('DOS (norm.)');
title('Edge vs bulk chain position'); legend('FontSize',8); grid on;

% Row 2 — LDOS maps and pair amplitude
subplot(3,4,5);
imagesc(1:N_chain, w, L2_chain);
axis xy; colorbar; colormap(gca,hot);
xlabel('Chain site'); ylabel('\omega');
title('Wire LDOS — proximity SC');
caxis([0, quantile(L2_chain(:),0.99)]);

subplot(3,4,6);
imagesc(1:N_chain, w, L3_chain);
axis xy; colorbar; colormap(gca,hot);
xlabel('Chain site'); ylabel('\omega');
title('Wire LDOS — topo (MZM at ends)');
caxis([0, quantile(L3_chain(:),0.99)]);

subplot(3,4,7);
imagesc(1:Nx, 1:Ny, F2_sub_2D);
axis xy equal tight; colorbar; colormap(gca,parula);
hold on;
plot([1,N_chain],[chain_row_bulk,chain_row_bulk],'w--','LineWidth',2);
xlabel('x'); ylabel('y');
title('|F(r)| substrate pair amplitude');
legend({'chain pos'},'FontSize',7,'Location','northeast');

subplot(3,4,8);
bar(1:N_chain, F2_wire, 'FaceColor',[0.2,0.5,0.9],'EdgeColor','none');
xlabel('Chain site i'); ylabel('|F(i)|');
title('Wire pair amplitude (proximity)'); grid on;

% Row 3 — Gamma scan, spectrum, geometry, MP
subplot(3,4,9);
yyaxis left;
plot(Gamma_vals, gap_vs_gamma, 'b-o','LineWidth',1.5,'MarkerSize',5);
ylabel('Quasiparticle gap');
yyaxis right;
plot(Gamma_vals, n_ingap_vs_gamma, 'r-s','LineWidth',1.5,'MarkerSize',5);
ylabel('# in-gap states');
xlabel('\Gamma'); title('Gap and in-gap states vs \Gamma'); grid on;

subplot(3,4,10);
E_all  = ev3.E;
tol3   = ex07_auto_tol(E_all);
in_gap = double(abs(E_all) < tol3);
scatter(1:numel(E_all), E_all, 4, in_gap, 'filled');
colormap(gca, [0.7 0.7 0.7; 1 0.2 0.2]);
caxis([0,1]); hold on;
if ~isempty(idx3)
    scatter(idx3, ev3.E(idx3), 30, 'r', 'filled');
end
yline(0,'k--','LineWidth',0.8);
xlabel('State index'); ylabel('E');
title(sprintf('Spectrum (%d in-gap, red)', numel(idx3)));
grid on;

subplot(3,4,11);
ax = gca; cla(ax);
% Substrate sites (scatter only — bonds via line objects for speed)
[XI,YI] = meshgrid(1:Nx, 1:Ny);
scatter(ax, XI(:), YI(:), 8, [0.8 0.8 0.9], 'filled'); hold(ax,'on');
% Vectorised bond drawing (2 line() calls instead of 1600 plot() calls)
xH = [reshape(XI(:,1:end-1),[],1), reshape(XI(:,2:end),[],1), nan(Ns2d-Ny,1)]';
yH = [reshape(YI(:,1:end-1),[],1), reshape(YI(:,2:end),[],1), nan(Ns2d-Ny,1)]';
xV = [reshape(XI(1:end-1,:),[],1), reshape(XI(2:end,:),[],1), nan(Ns2d-Nx,1)]';
yV = [reshape(YI(1:end-1,:),[],1), reshape(YI(2:end,:),[],1), nan(Ns2d-Nx,1)]';
line(ax, xH(:), yH(:), 'Color',[0.7 0.7 0.8], 'LineWidth',0.4);
line(ax, xV(:), yV(:), 'Color',[0.7 0.7 0.8], 'LineWidth',0.4);
% Wire sites and coupling
xc = 1:N_chain;
yc = chain_row_bulk * ones(1,N_chain);
xCoup = [xc; xc; nan(1,N_chain)];
yCoup = [yc; yc+0.5; nan(1,N_chain)];
line(ax, xCoup(:), yCoup(:), 'Color','g', 'LineWidth',1.2);
scatter(ax, xc, yc+0.5, 50, 'r', 'filled');
plot(ax, xc, yc+0.5, 'r-', 'LineWidth',1.5);
axis(ax,'equal'); axis(ax,[0.5 Nx+0.5 0.5 Ny+1]); axis(ax,'off');
title(ax, sprintf('Geometry: %dx%d substrate + %d-site wire', Nx, Ny, N_chain));
legend(ax,{'Substrate','','Wire',sprintf('\\Gamma=%.2f',Gamma)}, ...
    'Location','northeast','FontSize',7);

subplot(3,4,12);
if ~isempty(idx3) && any(mp_wire ~= 0)
    bar(1:N_chain, mp_wire, 'FaceColor',[0.2,0.6,0.3],'EdgeColor','none');
    hold on;
    xline(N_chain/2, 'r--', 'LineWidth',1.2);
    text(0.05, 0.92, sprintf('MP_{LR}=%.4f', MP_LR(1)), ...
        'Units','normalized','FontSize',8,'Color','r');
    lbl = 'topological MZM';
    if MP_LR(1) > 0, lbl = 'trivial / mixed'; end
    text(0.05, 0.80, lbl, 'Units','normalized','FontSize',7,'Color','r');
else
    text(0.15, 0.5, 'No in-gap states', 'Units','normalized','FontSize',9);
end
xlabel('Chain site'); ylabel('|MP(i)|');
title('Majorana polarization on wire'); grid on;

sgtitle(sprintf(['ex07 — 1D chain (N=%d) on 2D square lattice (%dx%d)\n' ...
    '\Gamma=%.2f   \Delta_s=%.2f   \lambda_R=%.2f   B_z=%.2f'], ...
    N_chain, Nx, Ny, Gamma, Delta_s, lambda_R, Bz), 'FontSize',11);

%% ═══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS  (must be at end of script in MATLAB)
%  ═══════════════════════════════════════════════════════════════════════════

function p = ex07_make_params(sub_sc, wire_topo, chain_row, Gamma, ...
                               Nx, Ny, N_chain, t_sub, t_wire, mu_sub, mu_wire, ...
                               Delta_s, lambda_R, Bz, eta)
    % Build UQMaT.Params for the chain1d_square2d heterostructure.
    %   sub_sc    : true → SC substrate (s-wave, Delta_s)
    %   wire_topo : true → Rashba+Zeeman wire (topological regime)

    if sub_sc
        sub_struct = struct('Nx',Nx,'Ny',Ny,'t',t_sub,'mu',mu_sub, ...
            'SC', struct('Delta_s',Delta_s,'sym','onsite'));
    else
        sub_struct = struct('Nx',Nx,'Ny',Ny,'t',t_sub,'mu',mu_sub);
    end

    if wire_topo
        wire_struct = struct('N',N_chain,'t',t_wire,'mu',mu_wire, ...
            'B',[0,0,Bz],'lambda_R',lambda_R);
    else
        wire_struct = struct('N',N_chain,'t',t_wire,'mu',mu_wire);
    end

    p = UQMaT.Params('system','chain1d_square2d', ...
        'hetero', struct( ...
            'substrate',      sub_struct, ...
            'chain',          wire_struct, ...
            'Gamma',          Gamma, ...
            'direction',      'x', ...
            'chain_row',      chain_row, ...
            'Gamma_disorder', 0.0), ...
        'eta', eta);
end

function tol = ex07_auto_tol(E)
    % Auto threshold for in-gap state detection.
    pos = E(E >  1e-10);
    neg = E(E < -1e-10);
    if ~isempty(pos) && ~isempty(neg)
        tol = min(pos(1), abs(neg(end))) * 0.5;
    else
        tol = 0.05 * (max(E) - min(E));
    end
end
