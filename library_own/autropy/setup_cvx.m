function setup_cvx

persistent cvx_setup_complete;

if (cvx_setup_complete)
    return;
end

global proteinroot;

str = [proteinroot,filesep,'pkg',filesep,'cvx',filesep,'cvx_setup.m'];

run(str);

cvx_setup_complete= 1;