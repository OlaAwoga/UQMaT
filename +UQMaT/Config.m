classdef Config
    % UQMaT.Config — Global configuration settings.
    %
    %   Controls MEX solver usage, verbosity, and other global options.
    %
    %   MEX solver control:
    %     UQMaT.Config.use_mex('auto')    % use MEX if compiled (default)
    %     UQMaT.Config.use_mex(true)      % always use MEX (error if not compiled)
    %     UQMaT.Config.use_mex(false)     % always use MATLAB fallback
    %     UQMaT.Config.mex_available()    % true if MEX files are present
    %     UQMaT.Config.mex_status()       % prints what is/isn't compiled
    %
    %   Verbosity:
    %     UQMaT.Config.verbose(true)      % print solver choices
    %     UQMaT.Config.verbose(false)     % silent (default)

    methods (Static)

        function use_mex(setting)
            % Set MEX usage: true | false | 'auto'
            if nargin == 0, setting = 'auto'; end
            p = UQMaT.Config.prefs_();
            if islogical(setting)
                p.mex = setting;
            elseif ischar(setting) || isstring(setting)
                switch lower(char(setting))
                    case 'auto',  p.mex = 'auto';
                    case 'on',    p.mex = true;
                    case 'off',   p.mex = false;
                    otherwise,    p.mex = 'auto';
                end
            end
            UQMaT.Config.set_prefs_(p);
        end

        function tf = mex_available(name)
            % Check if a specific MEX solver is available.
            %   UQMaT.Config.mex_available()           % any MEX
            %   UQMaT.Config.mex_available('eig')      % eigensolver
            %   UQMaT.Config.mex_available('green')    % GF solver
            %   UQMaT.Config.mex_available('kpm')      % KPM solver
            if nargin == 0
                tf = exist('uqmat_eig_mex','file') == 3 || ...
                     exist('uqmat_green_mex','file') == 3 || ...
                     exist('uqmat_kpm_mex','file') == 3;
            else
                tf = exist(['uqmat_' name '_mex'],'file') == 3;
            end
        end

        function tf = should_use_mex(name)
            % Internal: should we use MEX for this solver?
            p   = UQMaT.Config.prefs_();
            mex_setting = p.mex;
            if islogical(mex_setting)
                if mex_setting && ~UQMaT.Config.mex_available(name)
                    error('UQMaT:Config:mexNotCompiled', ...
                        ['MEX solver "%s" not compiled.\n' ...
                         'Run: UQMaT.Solver.Native.build_mex()'], name);
                end
                tf = mex_setting && UQMaT.Config.mex_available(name);
            else  % 'auto'
                tf = UQMaT.Config.mex_available(name);
            end
        end

        function mex_status()
            % Print which MEX solvers are compiled.
            fprintf('\n  UQMaT MEX solver status:\n');
            names = {'eig','green','kpm'};
            labels = {'Eigensolver (uqmat_eig_mex)', ...
                      'Green function (uqmat_green_mex)', ...
                      'KPM/Chebyshev (uqmat_kpm_mex)'};
            for k = 1:3
                tf = UQMaT.Config.mex_available(names{k});
                fprintf('  %s: %s\n', labels{k}, pick_(tf,'compiled ✓','not compiled'));
            end
            p = UQMaT.Config.prefs_();
            fprintf('  MEX setting: %s\n', num2str(p.mex));
            fprintf('  To compile: run UQMaT.Solver.Native.build_mex()\n\n');
        end

        function verbose(tf)
            % Set verbosity.
            p = UQMaT.Config.prefs_();
            p.verbose = tf;
            UQMaT.Config.set_prefs_(p);
        end

        function tf = is_verbose()
            p = UQMaT.Config.prefs_();
            tf = p.verbose;
        end

    end

    methods (Static, Access = private)
        function p = prefs_()
            persistent prefs
            if isempty(prefs)
                prefs.mex     = 'auto';
                prefs.verbose = false;
            end
            p = prefs;
        end

        function set_prefs_(p)
            % Store into persistent via a helper function
            UQMaT.Config.store_(p);
        end

        function store_(p)
            persistent prefs
            prefs = p;
        end
    end
end

function s = pick_(c,a,b)
    if c, s=a; else, s=b; end
end
