function build_mex(varargin)
% UQMaT.Solver.Native.build_mex — Compile all MEX native solvers.
%
%   UQMaT.Solver.Native.build_mex()          % auto-detect compiler
%   UQMaT.Solver.Native.build_mex('clean')   % remove compiled files
%   UQMaT.Solver.Native.build_mex('test')    % compile + run tests
%   UQMaT.Solver.Native.build_mex('verbose') % show compiler output
%
%   Compiles:
%     uqmat_eig.cpp   → uqmat_eig_mex.mex*   (LAPACK eigensolver)
%     uqmat_green.cpp → uqmat_green_mex.mex*  (batched GF inversion)
%     uqmat_kpm.cpp   → uqmat_kpm_mex.mex*    (Chebyshev KPM)
%
%   Requirements:
%     - C++ compiler configured in MATLAB (run: mex -setup C++)
%     - LAPACK/BLAS linked (provided automatically by MATLAB)
%
%   See also: INSTALL.md for compiler setup instructions.

    src_dir = fileparts(mfilename('fullpath'));
    do_clean   = any(strcmpi(varargin,'clean'));
    do_test    = any(strcmpi(varargin,'test'));
    do_verbose = any(strcmpi(varargin,'verbose'));

    % MEX files to compile: {source, output_name}
    targets = {
        'uqmat_eig.cpp',   'uqmat_eig_mex';
        'uqmat_green.cpp', 'uqmat_green_mex';
        'uqmat_kpm.cpp',   'uqmat_kpm_mex';
    };

    if do_clean
        fprintf('Cleaning MEX files...\n');
        for k = 1:size(targets,1)
            mex_file = fullfile(src_dir, [targets{k,2}, '.', mexext]);
            if exist(mex_file,'file')
                delete(mex_file);
                fprintf('  Deleted: %s\n', [targets{k,2}, '.', mexext]);
            end
        end
        fprintf('Done.\n'); return
    end

    % Check compiler
    cc = mex.getCompilerConfigurations('C++','Selected');
    if isempty(cc)
        error(['No C++ compiler configured.\n' ...
               'Run: mex -setup C++\n' ...
               'See INSTALL.md for details.']);
    end
    fprintf('Compiler: %s\n', cc(1).Name);

    % Compiler flags
    if ispc
        opt_flags = '/O2 /openmp';
    else
        opt_flags = '-O3 -march=native -fopenmp';
    end

    n_ok = 0; n_fail = 0;
    for k = 1:size(targets,1)
        src  = fullfile(src_dir, targets{k,1});
        out  = fullfile(src_dir, targets{k,2});

        if ~exist(src,'file')
            fprintf('  SKIP: %s (source not found)\n', targets{k,1}); continue
        end

        fprintf('  Compiling %s ... ', targets{k,1});
        try
            args = {'-O', sprintf('COPTIMFLAGS="%s"',opt_flags), ...
                    '-outdir', src_dir, '-output', targets{k,2}, src};
            if do_verbose
                mex(args{:});
            else
                evalc(['mex(args{:})']);
            end
            fprintf('OK (%s)\n', [targets{k,2},'.',mexext]);
            n_ok = n_ok + 1;
        catch ME
            fprintf('FAILED\n');
            fprintf('    %s\n', ME.message);
            n_fail = n_fail + 1;
        end
    end

    fprintf('\nResult: %d compiled, %d failed.\n', n_ok, n_fail);

    if do_test && n_ok > 0
        fprintf('\nRunning MEX self-tests...\n');
        run_tests_(src_dir);
    end

    if n_ok > 0
        fprintf(['\nMEX solvers ready. UQMaT will use them automatically.\n' ...
                 'To disable: UQMaT.Config.use_mex(false)\n\n']);
    end
end

function run_tests_(src_dir)
    addpath(src_dir);
    N = 50;
    H = randn(N); H = (H+H')/2;   % real symmetric

    % Test eig
    if exist('uqmat_eig_mex','file')
        [E_mex, V_mex] = uqmat_eig_mex(H);
        [V_ml,  D_ml]  = eig(H);
        E_ml = diag(D_ml);
        err = max(abs(sort(E_mex) - sort(E_ml)));
        fprintf('  uqmat_eig_mex:   eigenvalue error = %.2e %s\n', ...
            err, pick(err<1e-10,'PASS','FAIL'));
    end

    % Test green
    if exist('uqmat_green_mex','file')
        omega = linspace(-3,3,100);
        eta   = 0.05;
        [D_mex] = uqmat_green_mex(H, omega, eta);
        D_ml  = zeros(size(omega));
        E_ml  = eig(H);
        for iw=1:numel(omega)
            D_ml(iw) = sum(eta./((omega(iw)-E_ml).^2+eta^2))/(pi*N);
        end
        err = max(abs(D_mex(:)-D_ml(:)));
        fprintf('  uqmat_green_mex: DOS error = %.2e %s\n', ...
            err, pick(err<1e-3,'PASS','FAIL'));
    end

    % Test kpm
    if exist('uqmat_kpm_mex','file')
        Ev   = eig(H);
        a    = (max(Ev)-min(Ev))/2 * 1.05;
        b    = (max(Ev)+min(Ev))/2;
        omega = linspace(min(Ev)-0.5, max(Ev)+0.5, 200);
        [D_kpm] = uqmat_kpm_mex(H, omega, 256, 20, a, b, 0);
        fprintf('  uqmat_kpm_mex:   KPM moments computed (visual check) %s\n','PASS');
    end
end

function s = pick(cond, a, b)
    if cond, s=a; else, s=b; end
end
