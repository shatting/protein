clear classes

global proteinroot;
proteinroot = pwd;

global datadir;
datadir = [proteinroot,filesep,'data'];

% own libs
addpath([proteinroot,filesep,'classlib']);
addpath([proteinroot,filesep,'classlib',filesep,'utils']);
addpath([proteinroot,filesep,'library_own\hydros']);
addpath([proteinroot,filesep,'library_own\autropy']);
%addpath([proteinroot,filesep,'library_own\lga']);
addpath([proteinroot,filesep,'library_own\numhelices']);
addpath([proteinroot,filesep,'library_own\threefrags']);
%addpath([curdir,filesep,'main']);

% 3rd party libs
startup_lib('intlab');
%startup_lib('priority');
%startup_lib('DERIVESTsuite');
startup_lib('cvx');

disp('startup done.');