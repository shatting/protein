function createrunfile(p,probname,modprobname,algnum,maxit,data,clco,thfrags,newclasses,displaymode)

fid = fopen(probname,'w');

global remote_ampldir;
global remote_amploutdir;

 %[pathstr, filename, ext] = fileparts(probname);

n = length(data{p}.seq) - 3;

textv = ['## instructions for ampl \n'];
textv = [textv,'option solver "./knitroampl"; \n'];
textv = [textv,'option knitro_options "alg = ',num2str(algnum),' maxit = ',num2str(maxit),' maxcrossit = 5"; # maxcrossit=2 outlevel=1";\n'];
textv = [textv,' \n'];
textv = [textv,' model "',remote_ampldir,'/',modprobname,'_modfile.mod";\n'];
textv = [textv,' data "',remote_ampldir,'/',modprobname,'_protein.dat";\n'];
textv = [textv,' data "',remote_ampldir,'/_boundparams.dat";\n'];

if newclasses == 0
    textv = [textv,' data "',remote_ampldir,'/_potparams.dat";\n'];
else
    textv = [textv,' data "',remote_ampldir,'/_newpotparams.dat";\n'];
end

if thfrags == 1
    textv = [textv,' data "',remote_ampldir,'/_ft2.dat";\n'];
end
textv = [textv,' \n'];
textv = [textv,'solve; \n'];
textv = [textv,' \n'];
textv = [textv,' print {i in 1..',num2str(n+1),'}:c[i]  > ',remote_amploutdir,'/',modprobname,'_c.out;\n'];
textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];
%textv = [textv,' print {i in 1..',num2str(n),'}:ti[i]  > ',remote_amploutdir,'/',modprobname,'_ti.out;\n'];
%textv = [textv,' print {i in 1..',num2str(n),'}:tprime[i]  > ',remote_amploutdir,'/',modprobname,'_tprime.out;\n'];
%textv = [textv,' print {i in 1..',num2str(n+1),'}:s[i]  > ',remote_amploutdir,'/',modprobname,'_s.out;\n'];
%textv = [textv,' print {i in 1..',num2str(n),'}:phi[i]  > ',remote_amploutdir,'/',modprobname,'_phi.out;\n'];
if clco == 1
    textv = [textv,' print {i in 1..',num2str(n+2),'}:bond[i,1]  > ',remote_amploutdir,'/',modprobname,'_bond1.out;\n'];
    textv = [textv,' print {i in 1..',num2str(n+3),'}:calphas[i,1]  > ',remote_amploutdir,'/',modprobname,'_calphas1.out;\n'];
    textv = [textv,' print {i in 1..',num2str(n+2),'}:bond[i,2]  > ',remote_amploutdir,'/',modprobname,'_bond2.out;\n'];
    textv = [textv,' print {i in 1..',num2str(n+3),'}:calphas[i,2]  > ',remote_amploutdir,'/',modprobname,'_calphas2.out;\n'];
    textv = [textv,' print {i in 1..',num2str(n+2),'}:bond[i,3]  > ',remote_amploutdir,'/',modprobname,'_bond3.out;\n'];
    textv = [textv,' print {i in 1..',num2str(n+3),'}:calphas[i,3]  > ',remote_amploutdir,'/',modprobname,'_calphas3.out;\n'];
end

if thfrags == 1
    nfrags = 0;
    naa = size(data{p}.seq,1);
    for i = 1:(naa - 5 - 1)
        nfrags = nfrags + ((naa - 2) - (i + 3 + 1) + 1);
    end
    textv = [textv,' display thfragsmax;\n'];
    textv = [textv,' display sum{k in 1..',num2str(nfrags),'} wlnorm[k]*thfragsV[k]; \n'];
end

%textv = [textv,' print {i in 1..',num2str(n),'}:alph[i]  > ',remote_amploutdir,'/',modprobname,'_alph.out;\n'];
%textv = [textv,' print Vhy  > ',remote_amploutdir,'/',modprobname,'_Vhydro.out;\n'];
% textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];
% textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];
% textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];
% textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];
% textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];
% textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];
% textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];
% textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];
% textv = [textv,' print {i in 1..',num2str(n),'}:z[i]  > ',remote_amploutdir,'/',modprobname,'_z.out;\n'];

textv = [textv,' \n'];
% textv = [textv,'display calphasmean; \n'];
% textv = [textv,'display hphelper1; \n'];
% textv = [textv,'display hphelper2; \n'];
% textv = [textv,'display hp; \n'];
% textv = [textv,'display Vhydro; \n'];
% textv = [textv,'display cc; \n'];
% textv = [textv,' display ti;\n'];
% textv = [textv,'display tprime; \n'];
% textv = [textv,'display s; \n'];
% textv = [textv,' display c; \n'];
% textv = [textv,'display z; \n'];
% textv = [textv,'display phi; \n'];
%textv = [textv,'display z; \n'];
%textv = [textv,'display z; \n'];
% textv = [textv,'display firstccsum; \n'];
% textv = [textv,'  display fccsum;\n'];
if displaymode ~= 0
    if clco == 1
        textv = [textv,' display calphas;\n'];
        textv = [textv,' display c[1];\n'];
        textv = [textv,' display s[1];\n'];
        textv = [textv,'display bond; \n'];
        % textv = [textv,' display{i in 1..',num2str(n),',k in 1..',num2str(n),'}  (calphas[i,1]-calphas[k,1])^2 + (calphas[i,2]-calphas[k,2])^2 + (calphas[i,3]-calphas[k,3])^2;\n'];
         textv = [textv,'display {k in 1..',num2str(n),'} bond[k,2]*bond[k+1,3]-bond[k,3]*bond[k+1,2]; \n'];
        textv = [textv,'display a1; \n'];
        textv = [textv,'display {k in 1..',num2str(n),'} bond[k,3]*bond[k+1,1]-bond[k,1]*bond[k+1,3]; \n'];
        textv = [textv,'display a2; \n'];
        textv = [textv,'display {k in 1..',num2str(n),'} bond[k,1]*bond[k+1,2]-bond[k,2]*bond[k+1,1]; \n'];
        textv = [textv,' display a3;\n'];
        textv = [textv,'display deta;  \n'];
        textv = [textv,'display {k in 1..',num2str(n),'}  1/deta[k];  \n'];
        textv = [textv,'display detpartial; \n'];

        textv = [textv,'display {k in 1..',num2str(n),'} a1[k]*bond[k,2]*bond[k+1,3] + a2[k]*bond[k,3]*bond[k+1,1] + a3[k]*bond[k,1]*bond[k+1,2] \n'];
         textv = [textv,'       - a3[k]*bond[k,2]*bond[k+1,1] - a2[k]*bond[k,1]*bond[k+1,3] - a1[k]*bond[k,3]*bond[k+1,2]; \n'];

        if thfrags > 0
          textv = [textv,'display thfragsV; \n'];
          textv = [textv,'# display thfragsfeat;  \n'];   
           textv = [textv,'display{r in 1..',num2str(n+3),'}  sum{m in 1..3} bond[thfrags[r,1],m]*bond[thfrags[r,1]+1,m] - thfragsfeat1[r];      \n'];
           textv = [textv,'display thfragsfeat1;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  bond[thfrags[r,1],1]*tfr1[r] + bond[thfrags[r,1],2]*tfr2[r] + bond[thfrags[r,1],3]*tfr3[r] - thfragsfeat2[r];      \n'];
           textv = [textv,'display thfragsfeat2;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  bond[thfrags[r,1]+1,1]*tfr1[r] + bond[thfrags[r,1]+1,2]*tfr2[r] + bond[thfrags[r,1]+1,3]*tfr3[r] - thfragsfeat3[r];    \n'];
           textv = [textv,'display thfragsfeat3;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  sum{m in 1..3} bond[thfrags[r,1],m]*bond[thfrags[r,2],m] - thfragsfeat4[r];   \n'];
           textv = [textv,'display thfragsfeat4;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  sum{m in 1..3} bond[thfrags[r,1]+1,m]*bond[thfrags[r,2],m] - thfragsfeat5[r];  \n'];
           textv = [textv,'display thfragsfeat5;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  bond[thfrags[r,2],1]*tfr1[r] + bond[thfrags[r,2],2]*tfr2[r] + bond[thfrags[r,2],3]*tfr3[r] - thfragsfeat6[r];  \n'];
           textv = [textv,'display thfragsfeat6;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  sum{m in 1..3} bond[thfrags[r,1],m]*bond[thfrags[r,2]+1,m] - thfragsfeat7[r];   \n'];
           textv = [textv,'display thfragsfeat7;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  sum{m in 1..3} bond[thfrags[r,1]+1,m]*bond[thfrags[r,2]+1,m] - thfragsfeat8[r];   \n'];
           textv = [textv,'display thfragsfeat8;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  bond[thfrags[r,2]+1,1]*tfr1[r] + bond[thfrags[r,2]+1,2]*tfr2[r] + bond[thfrags[r,2]+1,3]*tfr3[r] - thfragsfeat9[r];   \n'];
           textv = [textv,'display thfragsfeat9;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  sum{m in 1..3} bond[thfrags[r,2],m]*bond[thfrags[r,2]+1,m] - thfragsfeat10[r];   \n'];
           textv = [textv,'display thfragsfeat10;  \n'];  
           textv = [textv,'display{r in 1..',num2str(n+3),'}  sqrt(tfr1[r]^2 + tfr2[r]^2 + tfr3[r]^2) - thfragsfeat11[r];    \n'];
           textv = [textv,'display thfragsfeat11;  \n'];  
        end
        textv = [textv,'  \n'];
        textv = [textv,'  \n'];
        % textv = [textv,' \n'];
        textv = [textv,'display c; \n'];
        textv = [textv,'display z; \n'];
        textv = [textv,'display cbuilder; \n'];
        textv = [textv,'display s; \n'];
        if newclasses == 0
        textv = [textv,'display {k in 1..',num2str(n),'} (2*cz[k]*z[k] + sz[k]*(1-z[k]^2))/(1 + z[k]^2)   ; \n'];
        end
        textv = [textv,'display ti; \n'];
        if newclasses == 0
        textv = [textv,'display {k in 1..',num2str(n),'} (- 2*sz[k]*z[k] + cz[k]*(1-z[k]^2))/(1 + z[k]^2); \n'];
        end
        textv = [textv,'display tprime; \n'];
        textv = [textv,'display s; \n'];
        textv = [textv,'display {k in 1..',num2str(n),'} bond[k + 2,1]^2 + bond[k + 2,2]^2 + bond[k + 2,3]^2; \n'];
        if newclasses == 1
            textv = [textv,'display t0; \n'];
            textv = [textv,'display tprime; \n'];
            textv = [textv,'display {k in 1..',num2str(n),'} z[k]/t0[k]; \n'];
            textv = [textv,'display ti; \n'];
            textv = [textv,'display {k in 1..',num2str(n),'} sqrt(1 - (z[k]/t0[k])^2); \n'];
        end

        textv = [textv,' \n'];
        %textv = [textv,'display cz; \n'];
        %textv = [textv,'display sz; \n'];

    end
end  %end if displaymode

fprintf(fid,textv);
fclose(fid);