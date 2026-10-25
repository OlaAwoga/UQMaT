function build_mex(varargin)
% build_mex — Compile UQMaT native MEX solvers.
%
%   build_mex()          % compile all (auto-detect compiler)
%   build_mex('clean')   % remove compiled files
%   build_mex('test')    % compile + run self-tests
%   build_mex('verbose') % show full compiler output
%
%   After compilation, UQMaT automatically uses the MEX solvers.
%   To control manually:
%     UQMaT.Config.use_mex('auto')    % MEX if compiled (default)
%     UQMaT.Config.use_mex(true)      % always MEX
%     UQMaT.Config.use_mex(false)     % always MATLAB
%     UQMaT.Config.mex_status()       % show what is compiled
%
%   See INSTALL.md for compiler requirements.

    % Add UQMaT to path if needed
    here = fileparts(mfilename('fullpath'));
    addpath(here);

    % Delegate to the Native package builder
    UQMaT.Solver.Native.build_mex(varargin{:});
end
