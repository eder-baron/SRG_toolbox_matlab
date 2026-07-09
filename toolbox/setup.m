%SETUP  Add the SRG Toolbox to the MATLAB search path.
%
%   Run this script once from the toolbox root directory before using the
%   toolbox, or add the call to your startup.m for permanent access.
%
%   Usage:
%       cd /path/to/srg_toolbox
%       run setup
%
%   To make the path permanent:
%       run setup
%       savepath
%
%   See also CONTENTS, SAVEPATH, PATHTOOL

root = fileparts(mfilename('fullpath'));
bg   = fullfile(root, 'bg');

if ~exist(bg, 'dir')
    error('setup:missingDir', ...
        'Directory bg/ not found. Run setup from the toolbox root.');
end

addpath(root);
addpath(bg);

fprintf('SRG Toolbox ready.\n');
fprintf('  Root : %s\n', root);
fprintf('  Utils: %s\n', bg);
fprintf('Type ''help Contents'' for a list of functions.\n');
