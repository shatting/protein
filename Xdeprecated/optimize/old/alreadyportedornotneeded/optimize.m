global proteinroot;
odir = pwd;
cd(proteinroot);

if (ispc)
    !optimize
else
    !./optimize.sh
end

cd(odir);