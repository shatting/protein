if ~exist([pwd,filesep,'startup_main.m'],'file'),
    disp('Warning: only do this in protein\matlab\ directory!');
    return;
end

disp('checking knitro(ktrlink)');
[x fval] = ktrlink(@(x)cos(x),1)

disp('startup: setting matlab path');

global proteinroot;
proteinroot = pwd;

% data dirs. we dont want them in path
%global rawdatadir; % RT* files
%rawdatadir = [proteinroot,filesep,'rawdata'];
global datadir; % data files with geometry info
datadir = [proteinroot,filesep,'database',filesep,'mats'];
global classdatadir; % 
classdatadir = [proteinroot,filesep,'main',filesep,'mats'];
global optdatadir; % 
optdatadir = [proteinroot,filesep,'optimize',filesep,'mats'];

global remote_proteinroot;
remote_proteinroot = '/users/shatting';

% ampl globals
global ampl_username;
ampl_username = 'shatting';
global local_ampldir;
local_ampldir = [proteinroot,filesep,'optimize',filesep,'ampl',filesep,'in'];
global remote_ampldir;
remote_ampldir = '/users/shatting/ampl/in';
global local_amploutdir;
local_amploutdir = [proteinroot,filesep,'optimize',filesep,'ampl',filesep,'out'];
global remote_amploutdir;
remote_amploutdir = '/users/shatting/ampl/out';


addpath(proteinroot);
addpath([proteinroot,filesep,'database']);
addpath([proteinroot,filesep,'geometry']);
addpath([proteinroot,filesep,'library']);
addpath([proteinroot,filesep,'main']);
addpath([proteinroot,filesep,'main',filesep,'interface']);
addpath([proteinroot,filesep,'optimize']);
addpath([proteinroot,filesep,'utils']);
addpath([proteinroot,filesep,'visual']);

cd([proteinroot,filesep,'library']);
startup_libs
cd(proteinroot)

disp_globals