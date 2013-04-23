function [ output_args ] = startup_lib( name )
%STARTUP_LIB 

persistent cvx_setup_complete;

global proteinroot;
libroot = [proteinroot,filesep,'library',filesep];
switch name
    case 'intlab'
        cd([libroot 'intlab']);
        startintlab;
        cd(proteinroot);
    case 'priority'
        addpath([libroot,'priority']);
    case 'DERIVESTsuite'
        addpath([libroot,'DERIVESTsuite']);
    case 'cvx'
        if (cvx_setup_complete)
            return;
        end
        run([proteinroot,filesep,'library',filesep,'cvx',filesep,'cvx_setup.m']);

        cvx_setup_complete= 1;

end

