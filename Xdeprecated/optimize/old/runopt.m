function  [GDT,LGArmsd,rmse] = runopt(data,class_info,newclasses,confus,probname,p,algnum,maxit,clco,hp,initials,avginitials,thfrags,zlu,displaymode,individualopt,s)
%rmse = runopt(data,class_info,newclasses,confus,probname,p,algnum,maxit,
%    clco,hp,initials,avginitials,thfrags,zlu,displaymode,individualopt,s)
% usually don't input s, unless want to try a specific secondary structure
% class sequence
%maxit: maximum number of iterations for ampl (default = 10000)
%alnum: options 1,2,or 3...
% algnum = 1: Knitro's Interior/Direct Algorithm
% algnum = 2: Knitro's Interior/CG Algorithm
% algnum = 3: Knitro's Active Set Algorithm
% alpha         max vgamma. if == 0 no upper bounds on Vgamma .
% beta          min vgamma. if == 0 no lower bounds on Vgamma.
% so, for alpha, choose a number like .95, which means that the value of
% Vgamma(i) must be less than the value that 95% of the sampled frequencies
% of class gamma 's Vgamma is less than
% createampl (not with createnewproblem)
% newconstraint: include/don't include the Vgamma(i) < Vgamma(j) constraint
% makefeasible : adjusts new constraint so that it is feasible
% clco: use close contacts constraint -- cc specifies amount allotted to be
% larger than actual cc. must be greater than 0 (cc = 0 means don't use
% this constraint)
% ex. cc = 1, 5, 10, 20, whatever
%try 137,323,194,191,190,188,175,169 -- all smallish
%if no class sequence provided, use actual sequence
% zlu (z lower,upper) = 1: lower and upper bounds for z
%                     = 0: only one joint bound for z (abs(z) < zboth)
% newclasses: 1, if hefg classes are used
% displaymode = 0: display less ampl output
% individualopt: used with boundchecker, otherwise set to 0. determines
% whether usual bounds or individual bounds on c and z are used 
% OUTPUT: 
% GDT: GDT_TS score from LGA
% LGArmsd: rmsd score from LGA
% rmse: rmsd score using register

if nargin < 17
   bond = rotatematrix(data{p}.bond); 
   if newclasses == 0
      s = seq2sclasses(bond); % protein s sequence
   else
       % skip, s is calculated with gamma, anyway
      %s = new_deal_hefg_interface(class_info,bond,data{p}.seq);
   end
   compare = 1; %if compare = 1: (default) runs getandcompare to compare 
                %new protein with old protein via rsmd, draw pics, and
                %make a pdb file
else
    compare = 0;
end

if hp == 1
    clco = 1;
end

if individualopt > 0
    zlu = 1;
end

if newclasses == 0
  gamma = seq2gammaclasses(data{p}.seq,class_info); % protein gamma sequence
else
  [s,gamma] = new_deal_hefg_interface(class_info,data{p}.bond,data{p}.seq);
end

global potparamswzupperlower;
if zlu == 1
    %if potparamswzupperlower = 1 %then we have the correct set of bounds
    %in potparams.dat
    if potparamswzupperlower ~= 1
        if newclasses == 0
        disp('need to make potparams again with zlower, zupper');
        disp('continuing with zboth');
        zlu = 0;
        end
    end
else % zlu = 0
    if potparamswzupperlower ~= 0
        disp('need to make potparams again with zboth');
        disp('continuing with zlower, zupper');
        zlu = 1;
    end
end
if newclasses == 1
    zlu = 1;
end

createampl(probname,data,p,algnum,maxit,gamma,s,clco,hp,initials,avginitials,zlu,thfrags,newclasses,displaymode,individualopt)
%    (probname, data, p, lookup, pot, alpha, beta, initials, avginitials, newconstraint,makefeasible, Vsort,algnum,maxit,bondls )
% s is the class sequence being sampled for p. default is the actual sequence

optimize
% this is going to annoy me:
%cd ginny;
rmse = 0;
GDT = 0;
LGArmsd = 0;

if compare == 1
    [GDT,LGArmsd,rmse] = getandcompare(p,data,probname,newclasses,confus,class_info,s);
end