function [geompred, X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = runknitro( priseq, sseq, gammaseq, potential, sbounds, sgbounds, infoneeded )
%RUNAMPL this is a rewrite of runopt.m
%
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

% if zlu == 1
%     %if potparamswzupperlower = 1 %then we have the correct set of bounds
%     %in potparams.dat
%     if potparamswzupperlower ~= 1
%         if newclasses == 0
%         disp('need to make potparams again with zlower, zupper');
%         disp('continuing with zboth');
%         zlu = 0;
%         end
%     end
% else % zlu = 0
%     if potparamswzupperlower ~= 0
%         disp('need to make potparams again with zboth');
%         disp('continuing with zlower, zupper');
%         zlu = 1;
%     end
% end
% if newclasses == 1
%     zlu = 1;
% end

nfrags = length(sseq);
dimx = nfrags + 1 + nfrags; % bond and torsion angles
ind_c = 1:nfrags+1;
ind_z = nfrags + 2 : dimx;

% generate objective
% x = [c_1, cp_1, z_1, ..., c_nfrags, cp_nfrags, z_nfrags]
% vis opticov.m: suff = suffstat(suff, double(feats), n, sysXto10([s
% gamma],[ns ng]));
potseq = sysXto10([sseq' gammaseq'],[infoneeded.pot_num_s infoneeded.pot_num_gamma]); % index sequence into potential variable
fobj = @(x) seq_pots(x,potseq,potential);

% generate bounds and starting values
% for bounds: intersect
% for X0: average c/cp mus from potential and use t mu
% TODO: X0 averageing could be done better as we have more information (cov) that we could use

% LB <= X <= UB
LB = zeros(dimx,1);
UB = LB;
X0 = UB;
% for c
for i = ind_c,    
    
    if (i>1 && i<=nfrags)
        frag_sgamma_i = potseq(i-1);
      	frag_sgamma_ip = potseq(i);
        % intersect bounds: c_i+1 \in \c(sgamma_i) and c_i \in \c(sgamma_i+1) 
        LB(i) = max(sgbounds.cplower(frag_sgamma_i),sgbounds.clower(frag_sgamma_ip));
        UB(i) = min(sgbounds.cpupper(frag_sgamma_i),sgbounds.cupper(frag_sgamma_ip));
       
        % take weighted average of means
        pot_i = potential(frag_sgamma_i);
        pot_ip = potential(frag_sgamma_ip);        
        X0(i) = (pot_i.mean(2)*pot_i.freq + pot_ip.mean(1)*pot_ip.freq) / (pot_ip.freq + pot_i.freq);        
    elseif i==1 % first
        frag_sgamma_i = potseq(1);        
        LB(i) = sgbounds.clower(frag_sgamma_i);
        UB(i) = sgbounds.cupper(frag_sgamma_i);
        X0(i) = potential(frag_sgamma_i).mean(1);
    else % last
        frag_sgamma_i = potseq(end);        
        LB(i) = sgbounds.cplower(frag_sgamma_i);
        UB(i) = sgbounds.cpupper(frag_sgamma_i);
        X0(i) = potential(frag_sgamma_i).mean(2);        
    end

    if (LB(i) > UB(i))
        dprintf('empty bound interval for bond angle c_%i',i)
    end
end

% bounds on z
% TODO: absolute bounds ("both")
for i = ind_z,    
    frag_sgamma_i = potseq(i-nfrags-1);
    LB(i) = sgbounds.zlower(frag_sgamma_i);
    UB(i) = sgbounds.zupper(frag_sgamma_i);   
    X0(i) = potential(frag_sgamma_i).mean(3);
end

% generate linear inequality constraints
A = [];
B = [];

% generate equality constraints
Aeq = [];
Beq = [];

%generate nonlinear (in)equality constraints
NONLINCON = [];
% TODO: elliptic constraints


%createampl(probname,data,p,algnum,maxit,gamma,s,clco,hp,initials,avginitials,zlu,thfrags,newclasses,displaymode,individualopt)
%    (probname, data, p, lookup, pot, alpha, beta, initials, avginitials, newconstraint,makefeasible, Vsort,algnum,maxit,bondls )
% s is the class sequence being sampled for p. default is the actual sequence

% KTRLINK(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS,knitroOptionsFile)
[X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = ktrlink(fobj,X0,A,B,Aeq,Beq,LB,UB,NONLINCON);

end
