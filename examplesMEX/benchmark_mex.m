%% benchmark_mex.m — MEX vs MATLAB timing benchmark for UQMaT solvers
%
%  Measures wall-clock time for eig(), Green's function, and KPM
%  at increasing system sizes. Plots speedup.
%
%  Requirements:
%    MEX compiled: run build_mex() first (needs C++ compiler).
%    See INSTALL.md for compiler setup.
%    Without MEX, the benchmark still runs MATLAB timings only.

root = mex_path();   % sets up path + MEX, returns UQMaT root directory

%% ── Check MEX availability ───────────────────────────────────────────────────
mex_eig   = UQMaT.Config.mex_available('eig');
mex_green = UQMaT.Config.mex_available('green');
mex_kpm   = UQMaT.Config.mex_available('kpm');
any_mex   = mex_eig || mex_green || mex_kpm;

fprintf('\n  MEX solver status:\n');
fprintf('  uqmat_eig_mex:   %s\n', pick(mex_eig,   'compiled ✓', 'NOT COMPILED'));
fprintf('  uqmat_green_mex: %s\n', pick(mex_green, 'compiled ✓', 'NOT COMPILED'));
fprintf('  uqmat_kpm_mex:   %s\n', pick(mex_kpm,   'compiled ✓', 'NOT COMPILED'));

if ~any_mex
    fprintf(['\n  No MEX solvers compiled.\n' ...
             '  Run build_mex() to compile (requires a C++ compiler).\n' ...
             '  See INSTALL.md for setup instructions.\n\n' ...
             '  Running MATLAB-only timings for reference...\n\n']);
end

%% ── System sizes ─────────────────────────────────────────────────────────────
N_vals = [50, 100, 200, 400];
if ~any_mex
    N_vals = [50, 100, 200];   % smaller range without MEX comparison
end
Nw     = 200;
omega  = linspace(-3, 3, Nw);
eta    = 0.03;
nrep   = 3;

t_eig_mex  = nan(1, numel(N_vals));
t_eig_ml   = nan(1, numel(N_vals));
t_grn_mex  = nan(1, numel(N_vals));
t_grn_ml   = nan(1, numel(N_vals));

fprintf('%-8s  %-14s  %-14s  %-14s  %-14s\n', ...
    'N', 'eig MEX (s)', 'eig MATLAB (s)', 'GF MEX (s)', 'GF MATLAB (s)');
fprintf('%s\n', repmat('-', 1, 70));

for in = 1:numel(N_vals)
    N   = N_vals(in);
    p   = UQMaT.Params('system','chain','N',N,'t',1,'mu',0);
    sys = UQMaT.System(p);
    H   = full(sys.H);

    % ── eig MEX ──────────────────────────────────────────────────────────────
    if mex_eig
        t = inf;
        for r = 1:nrep
            tic; uqmat_eig_mex(H); t = min(t, toc);
        end
        t_eig_mex(in) = t;
    end

    % ── eig MATLAB ───────────────────────────────────────────────────────────
    t = inf;
    for r = 1:nrep
        tic; eig(H); t = min(t, toc);
    end
    t_eig_ml(in) = t;

    % ── Green's function MEX ─────────────────────────────────────────────────
    if mex_green
        t = inf;
        for r = 1:nrep
            tic; uqmat_green_mex(H, omega, eta); t = min(t, toc);
        end
        t_grn_mex(in) = t;
    end

    % ── Green's function MATLAB ───────────────────────────────────────────────
    t = inf;
    for r = 1:nrep
        tic
        D_ = zeros(Nw, 1);
        for iw = 1:Nw
            z = omega(iw) + 1i*eta;
            G = (z*eye(N) - H) \ eye(N);
            D_(iw) = -imag(trace(G)) / (pi*N);
        end
        t = min(t, toc);
    end
    t_grn_ml(in) = t;

    % ── Print row ─────────────────────────────────────────────────────────────
    eig_mex_str = pick(mex_eig,   sprintf('%.4f', t_eig_mex(in)), 'N/A (compile MEX)');
    grn_mex_str = pick(mex_green, sprintf('%.4f', t_grn_mex(in)), 'N/A (compile MEX)');
    fprintf('%-8d  %-14s  %-14.4f  %-14s  %-14.4f\n', ...
        N, eig_mex_str, t_eig_ml(in), grn_mex_str, t_grn_ml(in));
end

%% ── Print speedup summary ────────────────────────────────────────────────────
fprintf('\n  Speedup summary:\n');
if mex_eig
    spd = t_eig_ml ./ t_eig_mex;
    fprintf('  eig:  '); fprintf('%.1fx  ', spd); fprintf('\n');
else
    fprintf('  eig:  MEX not compiled — run build_mex()\n');
end
if mex_green
    spd = t_grn_ml ./ t_grn_mex;
    fprintf('  GF:   '); fprintf('%.1fx  ', spd); fprintf('\n');
else
    fprintf('  GF:   MEX not compiled — run build_mex()\n');
end

%% ── Plot ─────────────────────────────────────────────────────────────────────
figure(1); clf;
set(gcf, 'Name', 'UQMaT MEX benchmark', 'Position', [100, 100, 1000, 380]);

subplot(1,3,1);
hold on;
if mex_eig
    loglog(N_vals, t_eig_mex, 'b-o', 'LineWidth',2, 'MarkerSize',8, ...
        'DisplayName','MEX (LAPACK)');
end
loglog(N_vals, t_eig_ml, 'r--s', 'LineWidth',2, 'MarkerSize',8, ...
    'DisplayName','MATLAB eig()');
xlabel('N'); ylabel('Time (s)');
title('Eigensolver scaling');
legend('FontSize',8, 'Location','northwest'); grid on;
if ~mex_eig
    text(0.05, 0.1, 'Compile MEX for comparison', ...
        'Units','normalized','FontSize',8,'Color','b');
end

subplot(1,3,2);
hold on;
if mex_green
    loglog(N_vals, t_grn_mex, 'b-o', 'LineWidth',2, 'MarkerSize',8, ...
        'DisplayName','MEX (LAPACK)');
end
loglog(N_vals, t_grn_ml, 'r--s', 'LineWidth',2, 'MarkerSize',8, ...
    'DisplayName','MATLAB loop');
xlabel('N'); ylabel('Time (s)');
title(sprintf('Green''s function (%d \\omega pts)', Nw));
legend('FontSize',8, 'Location','northwest'); grid on;
if ~mex_green
    text(0.05, 0.1, 'Compile MEX for comparison', ...
        'Units','normalized','FontSize',8,'Color','b');
end

subplot(1,3,3);
if any_mex
    hold on;
    if mex_eig
        spd_eig = t_eig_ml ./ t_eig_mex;
        plot(N_vals, spd_eig, 'b-o', 'LineWidth',2, 'MarkerSize',8, ...
            'DisplayName','eig speedup');
    end
    if mex_green
        spd_grn = t_grn_ml ./ t_grn_mex;
        plot(N_vals, spd_grn, 'r-s', 'LineWidth',2, 'MarkerSize',8, ...
            'DisplayName','GF speedup');
    end
    yline(1, 'k:', 'LineWidth',1);
    xlabel('N'); ylabel('Speedup (×)');
    title('MEX speedup over MATLAB');
    legend('FontSize',8); grid on;
else
    ax = gca; ax.Visible = 'off';
    text(0.1, 0.6, {'MEX solvers not compiled.','', ...
        'To enable:','  1. Install a C++ compiler','  2. Run build\_mex()',...
        '','See INSTALL.md'}, ...
        'Units','normalized','FontSize',10,'Color',[0.2 0.4 0.8]);
    title('Speedup (MEX not available)');
end

sgtitle('UQMaT MEX vs MATLAB benchmark', 'FontSize',12, 'FontWeight','bold');

fprintf('\nBenchmark complete.\n');
if ~any_mex
    fprintf(['To see MEX speedup:\n' ...
             '  1. Install C++ compiler (see INSTALL.md)\n' ...
             '  2. Run: build_mex()\n' ...
             '  3. Re-run: benchmark_mex\n\n']);
end

function s = pick(cond, a, b)
    if cond, s = a; else, s = b; end
end
