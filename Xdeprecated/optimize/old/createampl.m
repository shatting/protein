function createampl(probname, data, p, algnum, maxit,gamma,s,clco,hp,initials,avginitials,zlu,thfrags,newclasses,displaymode,individualopt)
% ex. createampl('newproblem',data,p,pot,class_info)
% prob name     problem name ex. 'newproblem' 
% s = the class sequence being sampled for p!!!
% data
% p             # of protein in data{}
% suffs         sufficient fragment statistics 
% pot           fragment potential from classification
% alpha         max vgamma. if == 0 no upper bounds.
% beta          min vgamma. if == 0 no lower bounds.
%               careful on those! if you change them, you have to set
%               redovgammas == 1
% initials      use actual initial values
%avginitials    starts opt. at avg. values of c, t, and tp for each class
% Vsort         sorted class potential values
% algnum     algnumber for ampl:
% algnum = 1: Knitro's Interior/Direct Algorithm
% algnum = 2: Knitro's Interior/CG Algorithm
% algnum = 3: Knitro's Active Set Algorithm
% maxit         maxiterations for ampl
%               = 0: (or just not there) no maximum - use knitro default
%               (10000, I think)
% createampl (not with createnewproblem)
% newconstraint: include/don't include the Vgamma(i) < Vgamma(j) constraint
% makefeasible : adjusts new constraint so that it is feasible
% bondls: use bond length constraints
% clco: use close contacts constraints
% clco = 0: don't use cc constraint
% clco = 1: use calphas constraint without relaxation 
% clco > 1 : allow sum of contacts to be greater than actual sum + cc
% so can't have sum + 1, or some plus anything between 0.0001 and 1
% inclusively - clco = 0,1, or anything larger than 1
% hp: include hydrophobic/philic constraints (from feat14)
% zlu (z lower,upper) = 1: lower and upper bounds for z
%                     = 0: only one joint bound for z (abs(z) < zboth)
% thfrags: use 3-frag constraint
% newclasses: 1, if hefg classes are used
% displaymode = 0: display less ampl output
% individualopt: used with boundchecker, otherwise set to 0. determines
% whether usual bounds or individual bounds on c and z are used 

global local_ampldir;

%fragclasses = classseq( data{p}.seq, lookup ); now s,gamma

%find value of close contacts constraint ( cc = 1/norm(xi - xk) forall xi, xk )
bond = rotatematrix(data{p}.bond); %rotates once (doesn't rotate bond2 yet) and normalizes
%bond = double(data{p}.bond);
%bond = bond./repmat(sqrt(sum(bond.^2,2)),1,3);
x = calphas(bond);
sizex = size(x,1);
actualcc = 0;
if hp > 0
    clcl = 1;
end
if clco > 0
 %actual cc is only sum of half, since symmetric
    for i = 1:sizex
        for k = 1:i-1
            actualcc = actualcc + 1/(.000001 + max(.0001,norm(x(i,:)-x(k,:),2))^2);
        end
    end
else
    actualcc = 0;
end

if thfrags == 1
    load thfragsstuff;
    pair2cl = tf.pair2cl;
    [thrfr,Vfrag] = get3fragpot(p,data,pair2cl,tf.potc);
    clear thfragsstff;
    n3f = 0;
    ninter = 1;
    for i = 1:(sizex - 5 - ninter)
       n3f = n3f + ((sizex - 2) - (i + 3 + ninter) + 1);
    end
    thrmax = max(max(pair2cl));
    if thrmax == 10 %this is sometimes = 10. don't know why. 
        return;
    end
else 
    thrfr = 0;
    n3f = 0;
    thrmax = 0;
    Vfrag = 0;
    pair2cl = 0;
end

%runopt(data,class_info,newclasses,confus,probname,p,algnum,maxit,
%    clco,hp,initials,avginitials,thfrags,zlu,displaymode,individualopt,s)
disp(['newclasses = ',num2str(newclasses)]);
disp(['p = ',num2str(p)]);
disp(['algnum = ',num2str(algnum)]);
disp(['maxit = ',num2str(maxit)]);
disp(['clco = ',num2str(clco)]);
disp(['hp = ',num2str(hp)]);
disp(['initials = ',num2str(initials)]);

disp(['avginitials = ',num2str(avginitials)]);
disp(['thfrags = ',num2str(thfrags)]);
disp(['zlu = ',num2str(zlu)]);
if individualopt == 1
disp('individualopt = just tighten z');
end
if individualopt == 2
disp('individualopt = just tighten c');
end
if individualopt == 3
disp('individualopt = tighten c and z');
end
    
disp('making protein.dat file');
createproteinfile(s,gamma,[local_ampldir,filesep,probname,'_protein.dat'],p,actualcc,hp,newclasses,thfrags,n3f,thrmax,thrfr,Vfrag,pair2cl);

% modfile
disp('making model file');
createmodfile([local_ampldir,filesep,probname,'_modfile.mod'], p, clco,hp,initials,avginitials,thfrags,n3f,zlu,newclasses,s,individualopt); 

% runfile
disp('making run file');
createrunfile(p,[local_ampldir,filesep,'test1','.run'],probname,algnum,maxit,data,clco,thfrags,newclasses,displaymode);
%1 at end just tells createrunfile that it was created with createampl
