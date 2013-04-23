% get globals, for if we were in a function
global geomdb;
global rawdatadir;
global datadir;

if (isempty(geomdb))
    if ask('load geomdb from file?')
        load([datadir,filesep,'geomdb.mat']);
        dprintf('global geomdb loaded.');    
    else              
        geomdb = data2geomdb([rawdatadir,filesep,'RT071127.mat']);
        if (ask('save geomdb?'))
            save([datadir,filesep,'geomdb.mat'],'geomdb');
            dprintf('global geomdb saved.');    
        end
    end
else
    dprintf('global geomdb already present, not loaded or generated.');
end
