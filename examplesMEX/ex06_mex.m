%% ex06_mex.m — MEX-accelerated version of examples/ex06_inhomogeneous_delta.m
%  Runs identical physics with MEX C++ solvers when compiled.
%  Falls back to MATLAB automatically if MEX not compiled.
%
%  Run build_mex() once before using this file.
%  See INSTALL.md for compiler setup.

root = mex_path();   % sets up path + MEX, returns UQMaT root directory

ex_file = fullfile(root, 'examples', 'ex06_inhomogeneous_delta.m');
if ~exist(ex_file, 'file')
    error('UQMaT:mex:notFound', ...
        'Cannot find: %s\nCheck that UQMaT is installed correctly.', ex_file);
end

%% ── Run WITH MEX ─────────────────────────────────────────────────────────────
fprintf('Running ex06_inhomogeneous_delta WITH MEX...\n');
t_mex = tic;

run(ex_file);

t_mex = toc(t_mex);
fprintf('[MEX]    %.3f s\n', t_mex);

%% ── Run WITHOUT MEX (MATLAB fallback) ────────────────────────────────────────
fprintf('\nRunning ex06_inhomogeneous_delta WITHOUT MEX (MATLAB fallback)...\n');
UQMaT.Config.use_mex(false);
t_ml = tic;

run(ex_file);

t_ml = toc(t_ml);
fprintf('[MATLAB] %.3f s\n', t_ml);

if UQMaT.Config.mex_available() && t_mex > 0
    fprintf('Speedup: %.1fx\n', t_ml / t_mex);
else
    fprintf('(MEX not compiled — no speedup comparison available)\n');
end

% Restore
UQMaT.Config.use_mex('auto');
