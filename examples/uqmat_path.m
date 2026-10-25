function root = uqmat_path()
% uqmat_path — Add UQMaT to MATLAB path.
%   Called at the top of every example script.
%   Returns the UQMaT root directory.
%
%   Note: Does NOT call clear classes (that wipes the return value).
%   To flush stale class definitions run: clear classes
%   in your MATLAB session before running examples, or use uqmat_setup.m.

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
             'Place examples/ inside the UQMaT root directory.'], here);
    end

    % Add to path — do NOT call clear classes here
    addpath(root);
    fprintf('UQMaT ready. Root: %s\n', root);
end
