%% uqmat_setup.m
%  Add UQMaT to the MATLAB path.
%  Call once per session: run('/full/path/to/UQMaT/uqmat_setup.m')
%
%  Or add permanently to startup.m:
%    edit(fullfile(userpath,'startup.m'))
%    % Add: run('/full/path/to/UQMaT/uqmat_setup.m')

root = fileparts(mfilename('fullpath'));  % folder containing this file
addpath(root);
clear classes   % flush stale class definitions (safe here — script, not function)
fprintf('UQMaT v1.0.0 ready. Root: %s\n', root);
