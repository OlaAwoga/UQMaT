function root = mex_path()
% mex_path — Set up UQMaT path and enable MEX solvers.
%   Called at the top of every examplesMEX script.
%   Returns the UQMaT root directory.
%
%   Note: Does NOT call clear classes (that wipes the return value).
%   If you need to flush stale class definitions, call clear classes
%   manually in your MATLAB session before running examples.

    here = fileparts(mfilename('fullpath'));

    % Walk up until +UQMaT is found
    root = here;
    for k = 1:5
        if exist(fullfile(root, '+UQMaT'), 'dir')
            break
        end
        root = fileparts(root);
    end

    if ~exist(fullfile(root, '+UQMaT'), 'dir')
        error('UQMaT:setup:notFound', ...
            ['Cannot find +UQMaT folder.\n' ...
             'Searched up from: %s\n' ...
             'Place examplesMEX/ inside the UQMaT root directory.'], here);
    end

    % Add to path — do NOT call clear classes here
    addpath(root);

    % MEX setup
    if UQMaT.Config.mex_available()
        UQMaT.Config.use_mex(true);
        UQMaT.Config.verbose(true);
        fprintf('[examplesMEX] MEX solvers enabled.\n');
    else
        warning('UQMaT:mex:notCompiled', ...
            ['MEX solvers not compiled. Run build_mex() first.\n' ...
             'Falling back to MATLAB implementation.\n' ...
             'See INSTALL.md for compiler setup.']);
        UQMaT.Config.use_mex(false);
    end

    fprintf('UQMaT root: %s\n', root);
end
