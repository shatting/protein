function createmodfile(probname, p,clco,hp,initials,avginitials,thfrags,n3f,zlu,newclasses,s,individualopt) %formatstring)
%createmodfile(problemname, p, initials, avginitials,bondls,cc)
%creates an ampl mod file for protein p ind RT*.mat
    
% zlu (z lower,upper) = 1: lower and upper bounds for z
%                     = 0: only one joint bound for z (abs(z) < zboth)
%if (nargin<3), formatstring = '%10.10f'; end

ninter = 1; %distance between 3frags

% can't use both types of initial values!
if initials == 1
    avginitials = 0;
end
if avginitials == 1
    initials = 0;
end

if thfrags > 0
    clco = 0;
end

%initials = 0;
%avginitials = 0;

%cc needs bonds!
if hp > 0
    clco = 1;
end
if clco > 0
    bondls = 1;
else
    bondls = 0;
end



%used for initial values: (if starting at actual values)
load RT071127.mat;
%bond = double(data{p}.bond);
%this is new -- rotate bond first (then bond -> geo -> bond works)
%bond = rot(double(bond));

%need a way to get geometry of a bond!!!
%geo = bond2geometry(bond); 
%c = [geo.c;geo.cp(end)];
%t= geo.t;
aaseq = data{p}.seq;
n = length(aaseq) - 3; %the number of 4-fragments

fid = fopen(probname,'w');
%[pathstr, filename, ext] = fileparts(probname);

textv = ['##Parameters!  \n'];
textv = [textv,' \n'];
textv = [textv,'# lengths and sequences for protein are stored in protein.dat \n'];
textv = [textv,'param nums;  ## number of secondary structure classes \n'];
textv = [textv,'param numg;  ## number of gamma classes  \n'];
textv = [textv,' \n'];
fprintf(fid,textv);
textv = [];
if (hp > 0 )  %|| thfrags > 0)
    textv = [textv,'param aaseq{i in 1..',num2str(n+3),'}; ## the amino acid sequence \n'];
    if thfrags > 0
        textv = [textv,'# unfortunately seems necessary for the thfrags for loop \n'];
       for i = 1:n+3
           textv = [textv,'let aaseq[',num2str(i),']:= ',num2str(aaseq(i)),'; \n'];
       end       
    end
    fprintf(fid,textv);
end

textv = [];
textv = [textv,'param seq{i in 1..',num2str(n),'}; ## the secondary structure class sequence \n'];
fprintf(fid,textv);
if newclasses == 1
    textz = [];
    textz = [textz,' # dont know why this is necessary    \n'];
        geo = bond2geo(data{p}.bond);
        seq = geo.tcl;
        for i = 1:size(geo.tcl,1)
            textz = [textz,'let seq[',num2str(i),'] := ',num2str(seq(i)),';   \n'];
        end
        fprintf(fid,textz);
end
textv = [];
textv = [textv,'param gseq{i in 1..',num2str(n),'}; ## the gamma class sequence \n'];
textv = [textv,'# these are like: seq[1] = 14, gseq[2] = 12,....all have fixed values between 1 and nums or numg \n'];
textv = [textv,' \n'];
textv = [textv,'# parameters in program.dat (always the same) \n'];
textv = [textv,'param ent{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,'param mu1{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,'param mu2{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,'param mu3{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,' \n'];
textv = [textv,'# I think R is supposed to be L^-1 or something, check! right now: R = L \n'];
textv = [textv,'param R11{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,'#param R12{s in 1..nums, g in 1..numg}; # got rid of, since just zeros\n'];
textv = [textv,'#param R13{s in 1..nums, g in 1..numg}; # got rid of, since just zeros\n'];
textv = [textv,'param R21{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,'param R22{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,'#param R23{s in 1..nums, g in 1..numg}; # got rid of, since just zeros\n'];
textv = [textv,'param R31{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,'param R32{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,'param R33{s in 1..nums, g in 1..numg}; \n'];
textv = [textv,' \n'];
textv = [textv,'# bound constraints for ci and ti, depend on si (s(i) = seq[i]) \n'];
if individualopt == 0
    textv = [textv,'param seqlowerc{i in 1..nums}; \n'];
    textv = [textv,'param sequpperc{i in 1..nums}; \n'];
    textv = [textv,'param seqlowercp{i in 1..nums}; \n'];
    textv = [textv,'param sequppercp{i in 1..nums}; \n'];
    if zlu == 1   
       textv = [textv,'param seqlowerz{i in 1..nums}; \n'];
       textv = [textv,'param sequpperz{i in 1..nums}; \n'];
    else
       textv = [textv,'param seqbothz{i in 1..nums}; \n'];
    end
elseif individualopt == 1
    textv = [textv,'param seqlowerc{i in 1..nums}; \n'];
    textv = [textv,'param sequpperc{i in 1..nums}; \n'];
    textv = [textv,'param seqlowercp{i in 1..nums}; \n'];
    textv = [textv,'param sequppercp{i in 1..nums}; \n'];
    if zlu == 1   
       textv = [textv,'param seqlowerz{i in 1..',num2str(n),'}; \n'];
       textv = [textv,'param sequpperz{i in 1..',num2str(n),'}; \n'];
    else
       textv = [textv,'param seqbothz{i in 1..',num2str(n),'}; \n'];
    end
elseif individualopt == 2
      textv = [textv,'param seqlowerc{i in 1..',num2str(n+1),'}; \n'];
    textv = [textv,'param sequpperc{i in 1..',num2str(n+1),'}; \n'];
    textv = [textv,'param seqlowercp{i in 1..',num2str(n+1),'}; \n'];
    textv = [textv,'param sequppercp{i in 1..',num2str(n+1),'}; \n'];
    if zlu == 1   
       textv = [textv,'param seqlowerz{i in 1..nums}; \n'];
       textv = [textv,'param sequpperz{i in 1..nums}; \n'];
    else
       textv = [textv,'param seqbothz{i in 1..nums}; \n'];
    end  
elseif individualopt == 3
    textv = [textv,'param seqlowerc{i in 1..',num2str(n+1),'}; \n'];
    textv = [textv,'param sequpperc{i in 1..',num2str(n+1),'}; \n'];
    textv = [textv,'param seqlowercp{i in 1..',num2str(n+1),'}; \n'];
    textv = [textv,'param sequppercp{i in 1..',num2str(n+1),'}; \n'];
    if zlu == 1   
       textv = [textv,'param seqlowerz{i in 1..',num2str(n),'}; \n'];
       textv = [textv,'param sequpperz{i in 1..',num2str(n),'}; \n'];
    else
       textv = [textv,'param seqbothz{i in 1..',num2str(n),'}; \n'];
    end
end
    
%textv = [textv,'param sequppert{i in 1..nums}; \n'];
textv = [textv,' \n'];
fprintf(fid,textv);
textv = [];
if clco > 0
   textv = ['#close contacts parameter: sum of close contacts constraint for target structure \n'];
   textv = [textv,'param cc; \n'];
   textv = [textv,'# alph, cz, sz depend on alpha and seq, and are used to find t,tprime from z \n'];
   textv = [textv,'#param alph{i in 1..',num2str(n),'}; \n'];
   if newclasses == 0
   textv = [textv,'param cz{i in 1..',num2str(n),'};   #cz_i is cos(alpha_i), alphas are midpoints of z ranges \n'];
   textv = [textv,'param sz{i in 1..',num2str(n),'};   #sz_i is sin(alpha_i), alphas are midpoints of z ranges \n'];
   end
   if hp > 0
       textv = [textv,'param hydroup; #upper bound for hydro constraint \n'];
       textv = [textv,'param hydrodown; #lower bound for hydro constraint \n'];
       textv = [textv,'param fitp; \n'];
       textv = [textv,'param hpL{i in 1..20}; \n'];
       textv = [textv,'param hpent{i in 1..20}; \n'];
       textv = [textv,'param hpmu{i in 1..20}; \n'];
   end
   fprintf(fid,textv);
   textv = [];
   
   %parameters for threefrags
   if thfrags > 0
       rmaxsq = 49;
       wlnormmax = 91125;
       textv = [textv,'#param n3f;  #number of 3-frags in this protein \n'];
       textv = [textv,'param n3fclass; # number of possible 3-frag classes \n'];
       textv = [textv,'param thfrags{i in 1..',num2str(n3f),', j in 1..2}; #all 3-frags in this protein \n'];
       textv = [textv,'param thfragsclass{i in 1..',num2str(n3f),'}; #class list of all 3-frags in this protein \n'];
       textv = [textv,'# param rmaxsq = 49; # have to figure out what the input of the weight function is\n'];
       textv = [textv,'# param wlnormmax = 91125; # have to figure out what the input of the weight function is\n'];
       textv = [textv,'param thfragsmax; #right now using the actual value as the maximal value of Vfrag \n'];
       textv = [textv,'# the parameters of the potential function Vfrag \n'];
       textv = [textv,'param thfragsmean{i in 1..n3fclass, j in 1..11}; \n'];
       textv = [textv,'param thfragsent{i in 1..n3fclass}; \n'];
       textv = [textv,'param thfragsL{i in 1..n3fclass, m in 1..11, n in 1..11}; \n'];
       textv = [textv,'#param pair2cl{i in 1..20,j in 1..20}; \n'];
       textv = [textv,'#seems to need to be defined in mod file bc of for loop \n'];
%        for i = 1:20
%            for j = 1:20
%                textv = [textv,'let pair2cl[',num2str(i),',',num2str(j),']:= ',num2str(pair2cl(i,j)),'; \n'];
%            end
%        end
   end % end if thfrags > 0
   
end % end if clco > 0
textv = [textv,' \n'];
textv = [textv,' \n'];
textv = [textv,'##Variables!!  \n'];
textv = [textv,'## (optimizing over x(i) = (c(i),c(i+1),z(i))) \n'];
textv = [textv,'var c{i in 1..',num2str(n+1),'}; \n'];
textv = [textv,'var z{i in 1..',num2str(n),'}; \n'];
if newclasses == 1
   textv = [textv,'## z is actually q, for the HEFG classes \n'];
end
fprintf(fid,textv);
textv = ['\n'];

global alphas;
if bondls == 1  
    %needed to find sz, cz (used to be done in protein.dat)
    if newclasses == 0
        alpha = [alphas(:,1:2),alphas(:,3:5)*pi/100];
        cz = [alpha(:,1:2),cos(alpha(:,3:5))];
        sz = [alpha(:,1:2),sin(alpha(:,3:5))];
        y = sys10toX(s,[3 3 3]);
        si = zeros(size(s));
        ci = zeros(size(s));
       % alphao = zeros(size(s));
        for i = 1:size(s,2)
           si(i)  = sz(y(i,2) + (y(i,1) - 1)*3,y(i,3) + 2);
           ci(i) =  cz(y(i,2) + (y(i,1) - 1)*3,y(i,3) + 2);
        %   alphao(i) = alpha(y(i,2) + (y(i,1) - 1)*3, y(i,3) + 2);
        end
    end
    
textv = [textv,' \n'];
fprintf(fid,textv);
 textz = ['\n'];
    textz = [textz,'var a1{k in 1..',num2str(n),'};\n'];
    textz = [textz,'var a2{k in 1..',num2str(n),'};\n'];
    textz = [textz,'var a3{k in 1..',num2str(n),'};\n'];
    textz = [textz,'var detpartial{k in 1..',num2str(n),'};\n'];
    textz = [textz,'var deta{k in 1..',num2str(n),'};\n'];
    textz = [textz,'var bond{k in 1..',num2str(n+2),',j in 1..3};\n'];  
    textz = [textz,'var calphas{i in 1..',num2str(n+3),', j in 1..3};\n']; 
%     textz = [textz,'#var tprime{i in 1..',num2str(n),'};\n'];
%     textz = [textz,'#var ti{i in 1..',num2str(n),'};\n'];
%     textz = [textz,'#var cbuilder{i in 1..',num2str(n+1),'};\n'];
%     textz = [textz,'#var s{i in 1..',num2str(n+1),'};\n'];
%     textz = [textz,'#var bu1{i in 1..',num2str(n),'};\n'];
%     textz = [textz,'#var bu2{i in 1..',num2str(n),'};\n'];
%     textz = [textz,'#var bu3{i in 1..',num2str(n),'};\n'];
    textz = [textz,'\n'];
     textz = [textz,'\n'];
     if newclasses == 0
            textz = [textz,'  ## cz and sz are cos, sin(alpha_i) \n'];
            textz = [textz,'  ## alpha_i is the midpoint of the range for phi_i \n'];
            textz = [textz,'  ## range for phi_i is dependent on s class  \n'];
            textz = [textz,'  ## so cz, sz are independent of the value of z, but \n'];
            textz = [textz,'  ## dependent on the s class - perfect! \n'];
            for i = 1:n
                textz = [textz,'  let cz[',num2str(i),'] := ',num2str(ci(i)),';\n'];
            end
            textz = [textz,'\n'];
            for i = 1:n
                textz = [textz,'  let sz[',num2str(i),'] := ',num2str(si(i)),';\n'];
            end
     end
    textz = [textz,'\n'];
    if newclasses == 0
        textz = [textz,'  var cbuilder{k in 1..',num2str(n+1),'} = c[k]/100;   ## need to define c8!!!    \n'];
    else 
        textz = [textz,'  var cbuilder{k in 1..',num2str(n+1),'} = c[k];  #/100;   ## need to define c8!!!    \n'];
    end
    textz = [textz,'  var s{k in 1..',num2str(n+1),'} = sqrt(1 - cbuilder[k]^2);  ## need to define s8!!!     \n'];
    
    if newclasses == 0
        textz = [textz,'  var tprime{k in 1..',num2str(n),'} = (- 2*sz[k]*z[k] + cz[k]*(1-z[k]^2))/(1 + z[k]^2); # t in geometryz2geometry        \n'];
        textz = [textz,'  var ti{k in 1..',num2str(n),'} = (2*cz[k]*z[k] + sz[k]*(1-z[k]^2))/(1 + z[k]^2);     \n'];
        textz = [textz,'     \n'];
    else
        textz = [textz,'  var t0{k in 1..',num2str(n),'} = s[k]*s[k+1];  \n'];
        textz = [textz,'  var tprime{k in 1..',num2str(n),'} =    \n'];
        textz = [textz,'        if seq[k] < 3 then (if seq[k] == 1 then sqrt(1 - (z[k]/max(t0[k],10^(-16)))^2)                             \n'];
        textz = [textz,'                                           else z[k]/max(t0[k],10^(-16))) \n'];
        textz = [textz,'                      else (if seq[k] == 3 then -sqrt(1 - (z[k]/max(t0[k],10^(-16)))^2)                         \n'];
        textz = [textz,'                                           else z[k]/max(t0[k],10^(-16)));  \n'];
        textz = [textz,'  var ti{k in 1..',num2str(n),'} =    \n'];
        textz = [textz,'        if seq[k] < 3 then (if seq[k] == 1 then z[k]/max(t0[k],10^(-16))                             \n'];
        textz = [textz,'                                           else sqrt(1 - (z[k]/max(t0[k],10^(-16)))^2))  \n'];
        textz = [textz,'                      else (if seq[k] == 3 then z[k]/max(t0[k],10^(-16))                        \n'];
        textz = [textz,'                                           else -sqrt(1 - (z[k]/max(t0[k],10^(-16)))^2));   \n'];
        textz = [textz,'     \n'];
        
    end
    
    textz = [textz,'  var bu1{k in 1..',num2str(n),'}= tprime[k]*s[k]*s[k+1];     \n'];
    textz = [textz,'  var bu2{k in 1..',num2str(n),'}=cbuilder[k]*cbuilder[k+1]+ti[k]*s[k]*s[k+1];      \n'];
    textz = [textz,'  var bu3{k in 1..',num2str(n),'}=cbuilder[k+1];     \n'];
    textz = [textz,'     \n'];
    textz = [textz,'  fix bond[1,1] :=  1;   \n'];
    textz = [textz,'  fix bond[1,2] :=  0;    \n'];
    textz = [textz,'  fix bond[1,3] :=  0;   \n'];
    textz = [textz,'  let bond[2,1] := cbuilder[1];    \n'];
    textz = [textz,'  let bond[2,2] := s[1];   \n'];
    textz = [textz,'  fix bond[2,3] := 0;   \n'];
    textz = [textz,'  fix calphas[1,1] :=  0;   \n'];
    textz = [textz,'  fix calphas[1,2] :=  0;   \n'];
    textz = [textz,'  fix calphas[1,3] :=  0;   \n'];
    textz = [textz,'  fix calphas[2,3] :=  0;   \n'];
    textz = [textz,'     \n'];
%     textz = [textz,'var flower{k in 1..',num2str(n),'};     \n'];
%     textz = [textz,'var kite{k in 1..',num2str(n),'};     \n'];
%     textz = [textz,'var star{k in 1..',num2str(n),'};     \n'];
%     textz = [textz,'var smiley{k in 1..',num2str(n),'};     \n'];
    fprintf(fid,textz);
    textz = [];
    

    textz = ['for {k in 1..',num2str(n),'} { \n'];
   textz = [textz,'   let a1[k] := bond[k,2]*bond[k+1,3]-bond[k,3]*bond[k+1,2]; \n'];
    textz = [textz,'   let a2[k] := bond[k,3]*bond[k+1,1]-bond[k,1]*bond[k+1,3]; \n'];
    textz = [textz,'   let a3[k] := bond[k,1]*bond[k+1,2]-bond[k,2]*bond[k+1,1]; \n'];
    textz = [textz,'\n'];
    textz = [textz,'                  \n'];
%     textz = [textz,'let flower[k] := (b3[k]/bond[k+1,3]) - (bond[k+1,1]/bond[k+1,3])*((b1[k]/a1[k]) - (a2[k]/a1[k])*((b2[k]/bond[k,2]) - (b1[k]*bond[k,1]/(a1[k]*bond[k,2])))/(1 - (a2[k]*bond[k,1]/((a1[k]*bond[k,2])))));                  \n'];
%     textz = [textz,'let kite[k] := (bond[k+1,2]/bond[k+1,3])*((b2[k]/bond[k,2]) - (b1[k]*bond[k,1]/(a1[k]*bond[k,2])))/(1 - (a2[k]*bond[k,1]/((a1[k]*bond[k,2]))));                  \n'];
%     textz = [textz,'let star[k] := ((a2[k]/a1[k])*((a3[k]*bond[k,1]/(a1[k]*bond[k,2])) - (bond[k+1,3]/bond[k,2])) - (1 - (a2[k]*bond[k,1]/(a1[k]*bond[k,2])))*a3[k]/a1[k])/(1 - (a2[k]*bond[k,1]/(a1[k]*bond[k,2])));                  \n'];
%     textz = [textz,'let smiley[k] := ((a3[k]*bond[k,1]/(a1[k]*bond[k,2])) - (bond[k,3]/bond[k,2]))/(1 - (a2[k]*bond[k,1])/(a1[k]*bond[k,2]));                  \n'];
%     textz = [textz,'let bond[k+2,3]:=  (flower[k] + kite[k])/(1 - (bond[k,1]/bond[k+1,3])*star[k] + (bond[k+1,2]*bond[k+1,3])*smiley[k]);                  \n'];
%     textz = [textz,'let bond[k+2,2]:= ((b2[k]/bond[k,2]) - (b1[k]*bond[k,1]/(a1[k]*bond[k,2])) + ((a3[k]*bond[k,1]/(a1[k]*bond[k,2])) - (bond[k,3]/bond[k,2]))*bond[k+2,3])/(1 - (a2[k]*bond[k,1]/(a1[k]*bond[k,2])));                  \n'];
%     textz = [textz,'let bond[k+2,1]:= (b1[k]/a1[k]) - (a2[k]/a1[k])*bond[k+2,2] - (a3[k]/a1[k])*bond[k+2,3];                  \n'];
%     textz = [textz,'  }; \n'];
    textz = [textz,'   let detpartial[k] := a1[k]*bond[k,2]*bond[k+1,3] + a2[k]*bond[k,3]*bond[k+1,1] + a3[k]*bond[k,1]*bond[k+1,2]\n'];
    textz = [textz,'           - a3[k]*bond[k,2]*bond[k+1,1] - a2[k]*bond[k,1]*bond[k+1,3] - a1[k]*bond[k,3]*bond[k+1,2]; \n'];
  %   textz = [textz,'   let deta[k] := if detpartial[k] <> 0 then 1/detpartial[k] else 10^(-14);   \n'];
     textz = [textz,'   let deta[k] := 1/max(detpartial[k],10^(-14)); \n'];
    
    textz = [textz,'   let bond[k+2,1] := (bu1[k]*bond[k,2]*bond[k+1,3] + a2[k]*bond[k,3]*bu3[k] + bu2[k]*bond[k+1,2]*a3[k] \n'];
    textz = [textz,'                                    - a3[k]*bond[k,2]*bu3[k] - a2[k]*bond[k+1,3]*bu2[k] - bond[k,3]*bond[k+1,2]*bu1[k])*deta[k];   \n'];
   
    textz = [textz,'   let bond[k+2,2] := (bu2[k]*a1[k]*bond[k+1,3] + bond[k+1,1]*bond[k,3]*bu1[k] + bu3[k]*bond[k,1]*a3[k]  \n'];
    textz = [textz,'           - a3[k]*bond[k+1,1]*bu2[k] - bond[k,3]*a1[k]*bu3[k] - bond[k,1]*bond[k+1,3]*bu1[k])*deta[k];     \n'];
   
    textz = [textz,'   let bond[k+2,3] :=  (bu3[k]*a1[k]*bond[k,2] + a2[k]*bond[k+1,1]*bu2[k] + bond[k+1,2]*bond[k,1]*bu1[k]   \n'];
    textz = [textz,'             - bond[k,2]*bond[k+1,1]*bu1[k] - bu2[k]*a1[k]*bond[k+1,2] - a2[k]*bond[k,1]*bu3[k])*deta[k];  \n'];
    %textz = [textz,'  }; \n']; 
    textz = [textz,'  # let calphas[k+1,1] := calphas[k,1] + bond[k,1];  \n'];
    textz = [textz,'  # let calphas[k+1,2] := calphas[k,2] + bond[k,2];  \n'];
    textz = [textz,'  # let calphas[k+1,3] := calphas[k,3] + bond[k,3];  \n'];
   textz = [textz,'                                                 \n'];
   textz = [textz,'};  \n'];
   textz = [textz,'subject to hi{k in 1..',num2str(n),'}: detpartial[k] == a1[k]*bond[k,2]*bond[k+1,3] + a2[k]*bond[k,3]*bond[k+1,1] + a3[k]*bond[k,1]*bond[k+1,2]     \n'];
   textz = [textz,'           - a3[k]*bond[k,2]*bond[k+1,1] - a2[k]*bond[k,1]*bond[k+1,3] - a1[k]*bond[k,3]*bond[k+1,2];    \n'];
   textz = [textz,'subject to hi2{k in 1..',num2str(n),'}: deta[k] == 1/max(detpartial[k],10^(-14));    \n'];
   textz = [textz,'subject to bb1{ k in 1..',num2str(n),'}: bond[k+2,1] == (bu1[k]*bond[k,2]*bond[k+1,3] + a2[k]*bond[k,3]*bu3[k] + bu2[k]*bond[k+1,2]*a3[k]      \n'];
   textz = [textz,'           - a3[k]*bond[k,2]*bu3[k] - a2[k]*bond[k+1,3]*bu2[k] - bond[k,3]*bond[k+1,2]*bu1[k])*deta[k];    \n'];
   textz = [textz,'subject to bb2{ k in 1..',num2str(n),'}: bond[k+2,2] == (bu2[k]*a1[k]*bond[k+1,3] + bond[k+1,1]*bond[k,3]*bu1[k] + bu3[k]*bond[k,1]*a3[k]      \n'];
   textz = [textz,'           - a3[k]*bond[k+1,1]*bu2[k] - bond[k,3]*a1[k]*bu3[k] - bond[k,1]*bond[k+1,3]*bu1[k])*deta[k];    \n'];
   textz = [textz,'subject to bb3{ k in 1..',num2str(n),'}: bond[k+2,3] ==  (bu3[k]*a1[k]*bond[k,2] + a2[k]*bond[k+1,1]*bu2[k] + bond[k+1,2]*bond[k,1]*bu1[k]     \n'];
   textz = [textz,'           - bond[k,2]*bond[k+1,1]*bu1[k] - bu2[k]*a1[k]*bond[k+1,2] - a2[k]*bond[k,1]*bu3[k])*deta[k];    \n'];
   
    fprintf(fid,textz);
    
end




textv = ['## need V(i) = norm(R[i](z(i) - mu(i)))^2 + ent(i) \n'];
textv = [textv,'# do the matrix multiplications and stuff termwise, since I have no idea how to get ampl to do them \n'];
textv = [textv,'var xminusmu1{i in 1..',num2str(n),'} = c[i]-mu1[seq[i],gseq[i]]; \n'];
textv = [textv,'var xminusmu2{i in 1..',num2str(n),'} = c[i+1]-mu2[seq[i],gseq[i]]; \n'];
textv = [textv,'var xminusmu3{i in 1..',num2str(n),'} = z[i]-mu3[seq[i],gseq[i]]; \n'];
textv = [textv,'var rtimesxminusmu1{i in 1..',num2str(n),'} = xminusmu1[i]*R11[seq[i],gseq[i]]; \n'];
textv = [textv,'var rtimesxminusmu2{i in 1..',num2str(n),'} = xminusmu1[i]*R21[seq[i],gseq[i]] \n'];
textv = [textv,'                                         + xminusmu2[i]*R22[seq[i],gseq[i]]; \n'];
textv = [textv,'var rtimesxminusmu3{i in 1..',num2str(n),'} = xminusmu1[i]*R31[seq[i],gseq[i]] \n'];
textv = [textv,'                                         + xminusmu2[i]*R32[seq[i],gseq[i]] \n'];
textv = [textv,'                                         + xminusmu3[i]*R33[seq[i],gseq[i]]; \n'];
textv = [textv,'##The potential variable \n'];
textv = [textv,'var Vsg{i in 1..',num2str(n),'} = rtimesxminusmu1[i]^2 + rtimesxminusmu2[i]^2 + rtimesxminusmu3[i]^2 + ent[seq[i],gseq[i]]; \n'];
textv = [textv,' \n'];
textv = [textv,'## bond length, c-alpha constraint variables could come here \n'];
fprintf(fid,textv);
textv = [];
if bondls == 1
    textz = ['# this is not correct yet. geometry2bond does not work yet - see thingstodo!! \n'];
    textz = [textz,'## bond length constraints!! (very tricky - have to be built up from the first and second bonds)\n'];
    textz = [textz,'\n'];
    textz = [textz,'  \n'];
    textz = [textz,'  \n'];
    %textz = [textz,'#var bn1{k in 1..',num2str(n),'};# = bn1b4det[k]/deta[k]   \n'];
    %textz = [textz,'#var bn2{k in 1..',num2str(n),'};# = bn2b4det[k]/deta[k];   \n'];
    %textz = [textz,'#var bn3{k in 1..',num2str(n),'};# = bn3b4det[k]/deta[k];  \n'];
    textz = [textz,'   \n'];
    textz = [textz,'#setting bnorm to 0 makes the program work, but I think it takes away the protection against\n'];
    textz = [textz,'#error propogation. but better to be a bit wrong than to be nothing at all...\n'];
    textz = [textz,'#the second and third columns of bond are a bit off without this, though \n'];
    textz = [textz,'#var bnorm{k in 1..',num2str(n),'} >=.0001;# sqrt(bn1[k]*bn1[k] + bn2[k]*bn2[k] + bn3[k]*bn3[k]);     \n'];
    textz = [textz,'     \n'];
    fprintf(fid,textz);
    textz = [];
end

if clco > 0
    texti = ['#only sum up half, bc symmetric\n'];
    texti = [texti,'var firstccsum{i in 1..',num2str(n+3),'} = sum{k in 1..i-1} 1/(.000001 + max(0.0001,(calphas[i,1]-calphas[k,1])^2 + (calphas[i,2]-calphas[k,2])^2 + (calphas[i,3]-calphas[k,3])^2));\n'];
    texti = [texti,'#var otherccsum{i in 1..',num2str(n+3),'} = sum{k in i+1..',num2str(n+3),'} 1/(.000001 + max(.0001,(calphas[i,1]-calphas[k,1])^2 + (calphas[i,2]-calphas[k,2])^2 + (calphas[i,3]-calphas[k,3])^2));\n'];
    texti = [texti,'var fccsum = sum{i in 1..',num2str(n+3),'} firstccsum[i];\n'];
    texti = [texti,'#var occsum = sum{i in 1..',num2str(n+3),'} otherccsum[i];\n'];
    texti = [texti,'     \n'];
    
    if hp > 0
        texti = [texti,'     \n'];
       texti = [texti,'# variables for hydrophobic/hydrophilic potential term: \n'];
       texti = [texti,'var calphasmean{j in 1..3} = sum{i in 1..',num2str(n+3),'} calphas[i,j]/',num2str(n+3),';  \n'];
       %texti = [texti,'var calphasmean2{j in 1..3} = sum{i in 1..',num2str(n+3),'} calphas[i,2];  \n'];
       %texti = [texti,'var calphasmean3{j in 1..3} = sum{i in 1..',num2str(n+3),'} calphas[i,3];  \n'];
       texti = [texti,'     \n'];
       texti = [texti,'var hphelper1{i in 1..',num2str(n+3),',j in 1..3} = calphas[i,j]/10 - calphasmean[j]/10; \n'];
       texti = [texti,'var hphelper2{i in 1..',num2str(n+3),'} = sum{j in 1..3} hphelper1[i,j]^2; \n'];
       texti = [texti,'var hp{i in 1..',num2str(n+3),'} = hphelper2[i]/fitp; \n'];
       texti = [texti,'     \n'];
       texti = [texti,'# the hydrophobic/hydrophilic potential term: \n'];
       texti = [texti,'# parameters are stored in protein.dat \n'];
       texti = [texti,'var Vhydro{i in 1..',num2str(n+3),'} = (hpL[aaseq[i]]*(hp[i] - hpmu[aaseq[i]]))^2 - hpent[aaseq[i]]/10;\n'];
       texti = [texti,'var Vhy = 100000*(sum{i in 1..',num2str(n+3),'} Vhydro[i]); \n'];
    end
    
    % 3-frag variables
    if thfrags > 0
       texti = [texti,' \n'];
       texti = [texti,'# Variables for 3-frag constraints! \n'];
       texti = [texti,'# var thfragsfeat{i in 1..',num2str(n3f),', j in 1..11}; \n'];
       texti = [texti,'var tfr1{i in 1..',num2str(n3f),'} = calphas[thfrags[i,2]+1,1] - calphas[thfrags[i,1]+1,1];  \n'];
       texti = [texti,'var tfr2{i in 1..',num2str(n3f),'} = calphas[thfrags[i,2]+1,2] - calphas[thfrags[i,1]+1,2];  \n'];
       texti = [texti,'var tfr3{i in 1..',num2str(n3f),'} = calphas[thfrags[i,2]+1,3] - calphas[thfrags[i,1]+1,3];  \n'];
       
%        texti = [texti,'# these seem to be needed for the loop to work \n'];
%        texti = [texti,'# thfrags(i,1) is index of first bond in first member of frag i \n'];
%        texti = [texti,'# thfrags(i,2) is index of first bond in second member of frag i \n'];
      % for i = 1:n3f
      %     texti = [texti,'let thfrags[',num2str(i),',1]:= ',num2str(thrfr(i,1)),';  \n'];
      %     texti = [texti,'let thfrags[',num2str(i),',2]:= ',num2str(thrfr(i,2)),';  \n'];
      % end 
       texti = [texti,'  \n'];
%        texti = [texti,'#let tfr1[1] := calphas[thfrags[1,2]+1,1] - calphas[thfrags[1,1]+1,1];    \n'];
%        texti = [texti,'#let tfr2[1] := calphas[thfrags[1,2]+1,2] - calphas[thfrags[1,1]+1,2];   \n'];
%        texti = [texti,'#let tfr3[1] := calphas[thfrags[1,2]+1,3] - calphas[thfrags[1,1]+1,3];      \n'];
       texti = [texti,'  \n'];
       texti = [texti,'  \n'];
       texti = [texti,'  \n'];
%        texti = [texti,'#let thfragsfeat[1,1]:= sum{k in 1..3} bond[thfrags[1,1],k]*bond[thfrags[1,1]+1,k];  \n'];
%        texti = [texti,'#let thfragsfeat[1,2]:= bond[thfrags[1,1],1]*tfr1[1] + bond[thfrags[1,1],2]*tfr2[1] + bond[thfrags[1,1],3]*tfr3[1]; \n'];
%        texti = [texti,'#let thfragsfeat[1,3]:= bond[thfrags[1,1]+1,1]*tfr1[1] + bond[thfrags[1,1]+1,2]*tfr2[1] + bond[thfrags[1,1]+1,3]*tfr3[1];  \n'];
%        texti = [texti,'#let thfragsfeat[1,4]:= sum{k in 1..3} bond[thfrags[1,1],k]*bond[thfrags[1,2],k];  \n'];
%        texti = [texti,'#let thfragsfeat[1,5]:= sum{k in 1..3} bond[thfrags[1,1]+1,k]*bond[thfrags[1,2],k];   \n'];
%        texti = [texti,'#let thfragsfeat[1,6]:= bond[thfrags[1,2],1]*tfr1[1] + bond[thfrags[1,2],2]*tfr2[1] + bond[thfrags[1,2],3]*tfr3[1];  \n'];
%        texti = [texti,'#let thfragsfeat[1,7]:= sum{k in 1..3} bond[thfrags[1,1],k]*bond[thfrags[1,2]+1,k];  \n'];
%        texti = [texti,'#let thfragsfeat[1,8]:= sum{k in 1..3} bond[thfrags[1,1]+1,k]*bond[thfrags[1,2]+1,k];  \n'];
%        texti = [texti,'#let thfragsfeat[1,9]:= bond[thfrags[1,2]+1,1]*tfr1[1] + bond[thfrags[1,2]+1,2]*tfr2[1] + bond[thfrags[1,2]+1,3]*tfr3[1];  \n'];
%        texti = [texti,'#let thfragsfeat[1,10]:= sum{k in 1..3} bond[thfrags[1,2],k]*bond[thfrags[1,2]+1,k];  \n'];
%        texti = [texti,'#let thfragsfeat[1,11]:= sqrt(tfr1[1]^2 + tfr2[1]^2 + tfr3[1]^2); \n'];
       texti = [texti,'  \n'];
       % texti = [texti,'param threefragschecki{k in 1..10};  \n'];
       % texti = [texti,'param threefragscheckj{k in 1..10};   \n'];
       % texti = [texti,'param r := 0;  \n'];
        texti = [texti,'var thfragsfeat1{r in 1..',num2str(n3f),'} = sum{m in 1..3} bond[thfrags[r,1],m]*bond[thfrags[r,1]+1,m];     \n'];
        texti = [texti,'var thfragsfeat2{r in 1..',num2str(n3f),'} = bond[thfrags[r,1],1]*tfr1[r] + bond[thfrags[r,1],2]*tfr2[r] + bond[thfrags[r,1],3]*tfr3[r];    \n'];
        texti = [texti,'var thfragsfeat3{r in 1..',num2str(n3f),'} = bond[thfrags[r,1]+1,1]*tfr1[r] + bond[thfrags[r,1]+1,2]*tfr2[r] + bond[thfrags[r,1]+1,3]*tfr3[r];     \n'];
        texti = [texti,'var thfragsfeat4{r in 1..',num2str(n3f),'} = sum{m in 1..3} bond[thfrags[r,1],m]*bond[thfrags[r,2],m];     \n'];
        texti = [texti,'var thfragsfeat5{r in 1..',num2str(n3f),'} = sum{m in 1..3} bond[thfrags[r,1]+1,m]*bond[thfrags[r,2],m];          \n'];
        texti = [texti,'var thfragsfeat6{r in 1..',num2str(n3f),'} = bond[thfrags[r,2],1]*tfr1[r] + bond[thfrags[r,2],2]*tfr2[r] + bond[thfrags[r,2],3]*tfr3[r];   \n'];
        texti = [texti,'var thfragsfeat7{r in 1..',num2str(n3f),'} = sum{m in 1..3} bond[thfrags[r,1],m]*bond[thfrags[r,2]+1,m];    \n'];
        texti = [texti,'var thfragsfeat8{r in 1..',num2str(n3f),'} = sum{m in 1..3} bond[thfrags[r,1]+1,m]*bond[thfrags[r,2]+1,m];   \n'];
        texti = [texti,'var thfragsfeat9{r in 1..',num2str(n3f),'} = bond[thfrags[r,2]+1,1]*tfr1[r] + bond[thfrags[r,2]+1,2]*tfr2[r] + bond[thfrags[r,2]+1,3]*tfr3[r];     \n'];
        texti = [texti,'var thfragsfeat10{r in 1..',num2str(n3f),'} = sum{m in 1..3} bond[thfrags[r,2],m]*bond[thfrags[r,2]+1,m];    \n'];
        texti = [texti,'var thfragsfeat11{r in 1..',num2str(n3f),'} = sqrt(tfr1[r]^2 + tfr2[r]^2 + tfr3[r]^2);    \n'];
        texti = [texti,'# param thfragsclass{r in 1..',num2str(n3f),'} = pair2cl[aaseq[thfrags[r,1]],aaseq[thfrags[r,2]]];    \n'];
%         texti = [texti,'subject to tf1{ r in 1..',num2str(n3f),'}: thfragsfeat1[r] == sum{m in 1..3} bond[thfrags[r,1],m]*bond[thfrags[r,1]+1,m]; \n'];
%         texti = [texti,'subject to tf2{ r in 1..',num2str(n3f),'}: thfragsfeat2[r] == bond[thfrags[r,1],1]*tfr1[thfrags[r,1]] + bond[thfrags[r,1],2]*tfr2[thfrags[r,1]] + bond[thfrags[r,1],3]*tfr3[thfrags[r,1]]; \n'];
%         texti = [texti,'subject to tf3{ r in 1..',num2str(n3f),'}: thfragsfeat3[r] == bond[thfrags[r,1]+1,1]*tfr1[thfrags[r,1]] + bond[thfrags[r,1]+1,2]*tfr2[thfrags[r,1]] + bond[thfrags[r,1]+1,3]*tfr3[thfrags[r,1]]; \n'];
%         texti = [texti,'subject to tf4{ r in 1..',num2str(n3f),'}: thfragsfeat4[r] == sum{m in 1..3} bond[thfrags[r,1],m]*bond[thfrags[r,2],m]; \n'];
%         texti = [texti,'subject to tf5{ r in 1..',num2str(n3f),'}: thfragsfeat5[r] == sum{m in 1..3} bond[thfrags[r,1]+1,m]*bond[thfrags[r,2],m];  \n'];
%         texti = [texti,'subject to tf6{ r in 1..',num2str(n3f),'}: thfragsfeat6[r] == bond[thfrags[r,2],1]*tfr1[thfrags[r,1]] + bond[thfrags[r,2],2]*tfr2[thfrags[r,1]] + bond[thfrags[r,2],3]*tfr3[thfrags[r,1]]; \n'];
%         texti = [texti,'subject to tf7{ r in 1..',num2str(n3f),'}: thfragsfeat7[r] == sum{m in 1..3} bond[thfrags[r,1],m]*bond[thfrags[r,2]+1,m];  \n'];
%         texti = [texti,'subject to tf8{ r in 1..',num2str(n3f),'}: thfragsfeat8[r] == sum{m in 1..3} bond[thfrags[r,1]+1,m]*bond[thfrags[r,2]+1,m]; \n'];
%         texti = [texti,'subject to tf9{ r in 1..',num2str(n3f),'}: thfragsfeat9[r] == bond[thfrags[r,2]+1,1]*tfr1[thfrags[r,1]] + bond[thfrags[r,2]+1,2]*tfr2[thfrags[r,1]] + bond[thfrags[r,2]+1,3]*tfr3[thfrags[r,1]]; \n'];
%         texti = [texti,'subject to tf10{ r in 1..',num2str(n3f),'}: thfragsfeat10[r] == sum{m in 1..3} bond[thfrags[r,2],m]*bond[thfrags[r,2]+1,m]; \n'];
%         texti = [texti,'subject to tf11{ r in 1..',num2str(n3f),'}: thfragsfeat11[r] == sqrt(tfr1[thfrags[r,1]]^2 + tfr2[thfrags[r,1]]^2 + tfr3[thfrags[r,1]]^2); \n'];
      
       %texti = [texti,'for {r in 1..',num2str(n3f),'} {  \n'];
       %texti = [texti,'for {i in 1..',num2str((n+3) - 5 - ninter),'} {  \n'];
       %texti = [texti,'for {j in i+',num2str(3 + ninter),'..',num2str(n+1),'} {  \n'];
           texti = [texti,'  \n'];
       %    texti = [texti,'let r := r+1;  \n'];
       %    texti = [texti,'let threefragschecki[r] := i;  \n'];
       %    texti = [texti,'let threefragscheckj[r] := j;  \n'];
           texti = [texti,'  \n'];
       %    texti = [texti,'let tfr1[r] := calphas[j+1,1] - calphas[i+1,1];    \n'];
       %    texti = [texti,'let tfr2[r] := calphas[j+1,2] - calphas[i+1,2];  \n'];
       %    texti = [texti,'let tfr3[r] := calphas[j+1,3] - calphas[i+1,3];   \n'];
           texti = [texti,'  \n'];
        %   texti = [texti,'let thfragsfeat[r,1]:= sum{m in 1..3} bond[i,m]*bond[i+1,m];   \n'];
        %   texti = [texti,'let thfragsfeat[r,2]:= bond[i,1]*tfr1[i] + bond[i,2]*tfr2[i] + bond[i,3]*tfr3[i]; \n'];
        %   texti = [texti,'let thfragsfeat[r,3]:= bond[i+1,1]*tfr1[i] + bond[i+1,2]*tfr2[i] + bond[i+1,3]*tfr3[i];  \n'];
        %   texti = [texti,'let thfragsfeat[r,4]:= sum{m in 1..3} bond[i,m]*bond[j,m];   \n'];
        %   texti = [texti,'let thfragsfeat[r,5]:= sum{m in 1..3} bond[i+1,m]*bond[j,m];   \n'];
        %   texti = [texti,'let thfragsfeat[r,6]:= bond[j,1]*tfr1[i] + bond[j,2]*tfr2[i] + bond[j,3]*tfr3[i];  \n'];
        %   texti = [texti,'let thfragsfeat[r,7]:= sum{m in 1..3} bond[i,m]*bond[j+1,m];  \n'];
        %   texti = [texti,'let thfragsfeat[r,8]:= sum{m in 1..3} bond[i+1,m]*bond[j+1,m]; \n'];
        %   texti = [texti,'let thfragsfeat[r,9]:= bond[j+1,1]*tfr1[i] + bond[j+1,2]*tfr2[i] + bond[j+1,3]*tfr3[i];  \n'];
        %   texti = [texti,'let thfragsfeat[r,10]:= sum{m in 1..3} bond[j,m]*bond[j+1,m];  \n'];
        %   texti = [texti,'let thfragsfeat[r,11]:= sqrt(tfr1[i]^2 + tfr2[i]^2 + tfr3[i]^2); \n'];
        %   texti = [texti,'let thfragsclass[r] := pair2cl[aaseq[i],aaseq[j]];  \n'];
      %  texti = [texti,' }; };  \n'];   
        texti = [texti,'  \n'];
        %texti = [texti,'var wl{i in 1..',num2str(n3f),'} = if ',num2str(rmaxsq),' >= thfragsfeat[i,11]^2 then (',num2str(rmaxsq),' - thfragsfeat[i,11]^2)^3 else 0;  \n'];
        texti = [texti,'var wl{i in 1..',num2str(n3f),'} = if 49 >= thfragsfeat11[i]^2 then (49 - thfragsfeat11[i]^2)^3 else 0;    \n'];
        texti = [texti,'var wlnorm{i in 1..',num2str(n3f),'} = (',num2str(wlnormmax),' - wl[i])/',num2str(wlnormmax),'; \n'];
        texti = [texti,'  \n'];
        texti = [texti,'var thfragsxminusmean1{i in 1..',num2str(n3f),'} = thfragsfeat1[i] - thfragsmean[thfragsclass[i],1];  \n'];
        texti = [texti,'var thfragsxminusmean2{i in 1..',num2str(n3f),'} = thfragsfeat2[i] - thfragsmean[thfragsclass[i],2];  \n'];
        texti = [texti,'var thfragsxminusmean3{i in 1..',num2str(n3f),'} = thfragsfeat3[i] - thfragsmean[thfragsclass[i],3];  \n'];
        texti = [texti,'var thfragsxminusmean4{i in 1..',num2str(n3f),'} = thfragsfeat4[i] - thfragsmean[thfragsclass[i],4];  \n'];
        texti = [texti,'var thfragsxminusmean5{i in 1..',num2str(n3f),'} = thfragsfeat5[i] - thfragsmean[thfragsclass[i],5];  \n'];
        texti = [texti,'var thfragsxminusmean6{i in 1..',num2str(n3f),'} = thfragsfeat6[i] - thfragsmean[thfragsclass[i],6];  \n'];
        texti = [texti,'var thfragsxminusmean7{i in 1..',num2str(n3f),'} = thfragsfeat7[i] - thfragsmean[thfragsclass[i],7];  \n'];
        texti = [texti,'var thfragsxminusmean8{i in 1..',num2str(n3f),'} = thfragsfeat8[i] - thfragsmean[thfragsclass[i],8];  \n'];
        texti = [texti,'var thfragsxminusmean9{i in 1..',num2str(n3f),'} = thfragsfeat9[i] - thfragsmean[thfragsclass[i],9];  \n'];
        texti = [texti,'var thfragsxminusmean10{i in 1..',num2str(n3f),'} = thfragsfeat10[i] - thfragsmean[thfragsclass[i],10];    \n'];
        texti = [texti,'var thfragsxminusmean11{i in 1..',num2str(n3f),'} = thfragsfeat11[i] - thfragsmean[thfragsclass[i],11];  \n'];
        texti = [texti,'  \n'];
        %texti = [texti,'var thfragsxminusmean{i in 1..',num2str(n3f),', j in 1..11} = thfragsfeat[i,j] - thfragsmean[thfragsclass[i],j];  \n'];
        %texti = [texti,'var thfragsLtimesmean{i in 1..',num2str(n3f),', j in 1..11} = sum{k in 1..11} thfragsL[thfragsclass[i],j,k]*thfragsxminusmean[i,k];  \n'];
        texti = [texti,'var thfragsLtimesmean1{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],1,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],1,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],1,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],1,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],1,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],1,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],1,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],1,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],1,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],1,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],1,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
        texti = [texti,'var thfragsLtimesmean2{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],2,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],2,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],2,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],2,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],2,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],2,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],2,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],2,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],2,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],2,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],2,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
       texti = [texti,'var thfragsLtimesmean3{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],3,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],3,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],3,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],3,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],3,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],3,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],3,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],3,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],3,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],3,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],3,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
       texti = [texti,'var thfragsLtimesmean4{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],4,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],4,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],4,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],4,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],4,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],4,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],4,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],4,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],4,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],4,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],4,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
       texti = [texti,'var thfragsLtimesmean5{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],5,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],5,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],5,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],5,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],5,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],5,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],5,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],5,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],5,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],5,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],5,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
       texti = [texti,'var thfragsLtimesmean6{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],6,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],6,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],6,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],6,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],6,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],6,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],6,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],6,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],6,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],6,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],6,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
       texti = [texti,'var thfragsLtimesmean7{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],7,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],7,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],7,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],7,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],7,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],7,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],7,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],7,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],7,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],7,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],7,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
       texti = [texti,'var thfragsLtimesmean8{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],8,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],8,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],8,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],8,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],8,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],8,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],8,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],8,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],8,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],8,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],8,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
       texti = [texti,'var thfragsLtimesmean9{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],9,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],9,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],9,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],9,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],9,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],9,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],9,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],9,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],9,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],9,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],9,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
       texti = [texti,'var thfragsLtimesmean10{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],10,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],10,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],10,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],10,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],10,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],10,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],10,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],10,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],10,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],10,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],10,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
       texti = [texti,'var thfragsLtimesmean11{i in 1..',num2str(n3f),'} = thfragsL[thfragsclass[i],11,1]*thfragsxminusmean1[i] + thfragsL[thfragsclass[i],11,2]*thfragsxminusmean2[i] + thfragsL[thfragsclass[i],11,3]*thfragsxminusmean3[i]  \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],11,4]*thfragsxminusmean4[i] + thfragsL[thfragsclass[i],11,5]*thfragsxminusmean5[i] + thfragsL[thfragsclass[i],11,6]*thfragsxminusmean6[i] + thfragsL[thfragsclass[i],11,7]*thfragsxminusmean7[i] \n'];
        texti = [texti,'                  + thfragsL[thfragsclass[i],11,8]*thfragsxminusmean8[i] + thfragsL[thfragsclass[i],11,9]*thfragsxminusmean9[i] + thfragsL[thfragsclass[i],11,10]*thfragsxminusmean10[i] + thfragsL[thfragsclass[i],11,11]*thfragsxminusmean11[i]; \n'];
        texti = [texti,'  \n'];
 
        texti = [texti,'var thfragsLtimesmeannormsquared{i in 1..',num2str(n3f),'} = thfragsLtimesmean1[i]^(2) + thfragsLtimesmean2[i]^(2) + thfragsLtimesmean3[i]^(2) + thfragsLtimesmean4[i]^(2)  \n'];
        texti = [texti,'                             + thfragsLtimesmean5[i]^(2) + thfragsLtimesmean6[i]^(2) + thfragsLtimesmean7[i]^(2) + thfragsLtimesmean8[i]^(2)  \n'];
        texti = [texti,'                             + thfragsLtimesmean9[i]^(2) + thfragsLtimesmean10[i]^(2) + thfragsLtimesmean11[i]^(2);  \n'];
       % texti = [texti,'var thfragsV{i in 1..',num2str(n3f),'} = thfragsent[thfragsclass[i]] + sum{k in 1..11} thfragsLtimesmean[i,k]^2;  \n'];
        texti = [texti,'var thfragsV{i in 1..',num2str(n3f),'} = thfragsent[thfragsclass[i]] + thfragsLtimesmeannormsquared[i];  \n'];
        texti = [texti,'  \n'];
        texti = [texti,'  \n'];
        texti = [texti,'  \n'];
 
    end % end if thfrags > 0
        
    fprintf(fid,texti);   
    texti = [];
end

if ((newclasses == 1) && (clco == 1))
    textu = ['# to prevent a whole columnn of 0s in bond matrix \n'];
    textu = [textu,'var bond3sum = sum{k in 1..',num2str(n+1),'} bond[k,3]; \n'];
    fprintf(fid,textu);    
    textu = [];
end

textv = [textv,' \n'];
textv = [textv,'##Problem! \n'];
textv = [textv,' \n'];
textv = [textv,'# minimize the sum of all the Vsg:  \n'];
textv = [textv,'minimize Problem: sum{i in 1..',num2str(n),'} Vsg[i]; \n'];

textv = [textv,' \n'];
textv = [textv,'##Constraints... \n'];
textv = [textv,'#these are based on heuristic bounds (using findsbounds), not the defined bounds. \n'];
textv = [textv,' \n'];
if ((individualopt ~= 2) && (individualopt ~= 3)) %if normal c bounds
    textv = [textv,'subject to cbounds{i in 1..',num2str(n),'}: seqlowerc[seq[i]] <= c[i] <= sequpperc[seq[i]]; \n'];
    textv = [textv,'subject to cpbounds{i in 2..',num2str(n+1),'}: seqlowercp[seq[i-1]] <= c[i] <= sequppercp[seq[i-1]]; \n'];
else %individual c bounds
    textv = [textv,'subject to cbounds{i in 1..',num2str(n),'}: seqlowerc[i] <= c[i] <= sequpperc[i]; \n'];
    textv = [textv,'subject to cpbounds{i in 2..',num2str(n+1),'}: seqlowercp[i-1] <= c[i] <= sequppercp[i-1]; \n'];
end
    
    
if zlu == 1
    if ((individualopt ~= 1) && (individualopt ~= 3)) %if normal z bounds
       textv = [textv,'subject to zbounds{i in 1..',num2str(n),'}: seqlowerz[seq[i]] <= z[i] <= sequpperz[seq[i]]; \n'];
    else
        textv = [textv,'subject to zbounds{i in 1..',num2str(n),'}: seqlowerz[i] <= z[i] <= sequpperz[i]; \n'];
    end      
else
    if ((individualopt ~= 1) && (individualopt ~= 3)) %if normal z bounds
        textv = [textv,'subject to zbounds{i in 1..',num2str(n),'}: -seqbothz[seq[i]] <= z[i] <= seqbothz[seq[i]]; \n'];
    else
        textv = [textv,'subject to zbounds{i in 1..',num2str(n),'}: -seqbothz[i] <= z[i] <= seqbothz[i]; \n'];
    end  
end
textv = [textv,' \n'];
if hp == 1
textv = [textv,'subject to hydroconstraint: hydrodown <= Vhy <= hydroup; \n'];
textv = [textv,' \n'];
end
if clco > 0
    textv = [textv,'#more constraints for c-alphas,... come here\n'];
    textv = [textv,'subject to aa1{k in 1..',num2str(n),'}: a1[k] == bond[k,2]*bond[k+1,3]-bond[k,3]*bond[k+1,2];  \n'];
    textv = [textv,'subject to aa2{k in 1..',num2str(n),'}: a2[k] == bond[k,3]*bond[k+1,1]-bond[k,1]*bond[k+1,3]; \n'];
    textv = [textv,'subject to aa3{k in 1..',num2str(n),'}: a3[k] == bond[k,1]*bond[k+1,2]-bond[k,2]*bond[k+1,1];     \n'];
    textv = [textv,' \n'];
   % if thfrags == 1
   % textv = [textv,'subject to detprob{k in 1..',num2str(n),'}: detpartial[k] == a1[k]*bond[k,2]*bond[k+1,3] + a2[k]*bond[k,3]*bond[k+1,1] + a3[k]*bond[k,1]*bond[k+1,2] \n'];
   % textv = [textv,'           - a3[k]*bond[k,2]*bond[k+1,1] - a2[k]*bond[k,1]*bond[k+1,3] - a1[k]*bond[k,3]*bond[k+1,2];   \n'];
   % textv = [textv,'subject to hatedet{k in 1..',num2str(n),'}: deta[k] == if detpartial[k] <> 0 then 1/detpartial[k] else 10^(-14);  \n'];
   % textv = [textv,' \n'];
   % end
    textv = [textv,'subject to bond1: bond[2,1] == cbuilder[1]; \n'];
    textv = [textv,'subject to bond2: bond[2,2] == s[1]; \n'];
    textv = [textv,'subject to calpac{k in 1..',num2str(n+2),',j in 1..3}: calphas[k+1,j] == bond[k,j] + calphas[k,j];  \n'];
    if newclasses == 1
        textv = [textv,'# to prevent a whole column of 0s in the bond matrix \n'];
        textv = [textv,'subject to bondextra1: abs(bond3sum) >= .001; \n'];
        %textv = [textv,'subject to bondextra2: bond3sum <= -.001; \n'];
    end
end

fprintf(fid,textv);
textv = [];
if bondls == 1
    textz = [];
%     textz = [textz,'##bond length constraints\n'];
%     textz = [textz,'subject to bond11 : bond[1,1] == 1;\n'];
%     textz = [textz,'subject to bond12 : bond[1,2] == 0;\n'];
%     textz = [textz,'subject to bond13 : bond[1,3] == 0;\n'];
%     textz = [textz,'subject to bond21 : bond[2,1] == cbuilder[1];\n'];
%     textz = [textz,'subject to bond22 : bond[2,2] == s[1];\n'];
%     textz = [textz,'subject to bond23 : bond[2,3] == 0;\n'];
    textz = [textz,'    \n'];
    textz = [textz,'    \n'];
    textz = [textz,'#defines the bond lengths   \n'];
     textz = [textz,'    \n'];
%    textz = [textz,'subject to detaconstraint{k in 1..7}: deta[k] == a1[k]*bond[k,2]*bond[k+1,3] + a2[k]*bond[k,3]*bond[k+1,1] + a3[k]*bond[k,1]*bond[k+1,2]    \n'];
%    textz = [textz,'      - a3[k]*bond[k,2]*bond[k+1,1] - a2[k]*bond[k,1]*bond[k+1,3] - a1[k]*bond[k,3]*bond[k+1,2];     \n'];
%    textz = [textz,'subject to detnotzero{i in 1..7}:     \n'];
%    textz = [textz,'       if deta[i] == 0\n'];
%    textz = [textz,'           then deta[i] = 10^(-14); \n'];
%    textz = [textz,'subject to bondnotzero{i in 3..9,j in 1..3}:    \n'];
%    textz = [textz,'   if bond[i,j] == 0  \n'];
%    textz = [textz,'       then bond[i,j] = 10^(-14); \n'];
%    textz = [textz,'    \n'];
%    textz = [textz,'subject to bd1{ k in 1..',num2str(n),'}: bn1[k] == bn1b4det[k]/deta[k];    \n'];
%    textz = [textz,'subject to bd2{ k in 1..',num2str(n),'}: bn2[k] == bn2b4det[k]/deta[k];    \n'];
%    textz = [textz,'subject to bd3{ k in 1..',num2str(n),'}: bn3[k] == bn3b4det[k]/deta[k];    \n'];
%    textz = [textz,'    \n'];
    textz = [textz,'    \n'];
    
    
    textz = [textz,'#restrictions on bond lengths: (so I can actually get rid of definitions of bonds)    \n'];
    textz = [textz,'subject to lowerbondbound{ i in 1..',num2str(n),'}: bond[i + 2,1]^2 + bond[i + 2,2]^2 + bond[i + 2,3]^2  >= .95^2;     \n'];
    textz = [textz,'subject to upperbondbound{ i in 1..',num2str(n),'}: bond[i + 2,1]^2 + bond[i + 2,2]^2 + bond[i + 2,3]^2  <= 1.05^2;     \n'];
    textz = [textz,'# I think norm(bond(i+2)) is always 1, by the way geometry2bond was defined.     \n'];
    textz = [textz,'    \n'];
    fprintf(fid,textz);
    textz = [];
end
if clco > 0
    textz = ['# subject to calpa1{j in 1..3}: calphas[1,j] == 0;\n'];
    textz = [textz,'#close contact constraints:\n'];
    textz = [textz,'#maybe it would be smarter to only sum of the first half of cc, too-- would\n'];
    textz = [textz,'#probably save time in the long run -- okay, that is what i am doing now - thanks for the tip, me! \n'];
    if clco > 1
       textz = [textz,'subject to ccupper: fccsum  <= cc + ',num2str(clco),';\n']; %   sum{i in 1..',num2str(n),'} \n'];
       %textz = [textz,'           sum{k in 1..',num2str(n),'} ((calphas[i,1]-calphas[k,1])^2 + (calphas[i,2]-calphas[k,2])^2 + (calphas[i,3]-calphas[k,3])^2)^(-1/2) <= cc + 500;\n'];
    else
        textz = [textz,'subject to ccupper: fccsum  <= cc;\n']; %   sum{i in 1..',num2str(n),'} \n'];
    end
    fprintf(fid,textz);
    textz = [];
end

if thfrags == 1
    textz = ['# 3-frag constraints! \n'];
    textz = [textz,'# need to add weight to sum!!! \n'];
    textz = [textz,'subject to thfragsconstraint: sum{k in 1..',num2str(n3f),'} wlnorm[k]*thfragsV[k] <= thfragsmax + 10;  \n'];
    textz = [textz,' \n'];
    fprintf(fid,textz);
    textz = [];
end



textv = [textv,' \n'];
textv = [textv,'# initial values come here! \n'];
fprintf(fid,textv);
textv = [];

global alphas;
% bond initials
   
if (initials)
     %% set the cs to start very close to their original values
     if newclasses == 0
         geo = bond2geometry(data{p}.bond);
         geomz = geometry2geometryz(geo, alphas);
         c = 100*geomz.co;
         z = geomz.zo;
     else % newclasses == 1
         
         c = [geo.co;geo.cpo(end)];
         z = geo.tqu;
     end
     
     textv = ['# initial values - set very close to actual values.\n'];
     for i = 1:size(c,1)
         textv = [textv,'let c[',num2str(i),']:=',num2str(c(i)),';\n'];
     end
     fprintf(fid,textv);
     
     %% set the zs to start very close to their original values
     textv = [' \n'];
     for i = 1:size(z,1)
         textv = [textv,'let z[',num2str(i),']:=',num2str(z(i)),';\n'];
     end
    fprintf(fid,textv);
    textv = [];
end 

if (avginitials)
    load classavg;
         % these are the average values of c, cp, and z for each class
         if newclasses == 1
             avgc = classavgs(:,1);
             avgcp = classavgs(:,2);
             avgz = classavgs(:,3);
         else
             avgc = 100*classavgs(:,1);
             avgcp = 100*classavgs(:,2);
             avgz = classavgs(:,3);
         end
         textv = ['# using avg value of c, cp, and z for each class \n'];
         textv = [textv,'# not really sure what to do about c, cp, so just averaging them \n'];
         numclasses = size(avgc,1); % size(avgc) = num classes!
         avgccp = zeros(numclasses,1);
         for i = 2:numclasses-1 
             avgccp(i) = (avgc(i) + avgcp(i-1))/2;
         end
         avgccp(1) = avgc(1);
         avgccp(numclasses) = avgcp(numclasses);
         if newclasses == 1
             seq = geo.tcl;
         else 
             seq = seq2sclasses(data{p}.bond);
         end   
         for i = 1:length(seq)
             s = seq(i);
             textv = [textv,'let c[',num2str(i),'] := ',num2str(avgccp(s)),'; \n'];
             textv = [textv,'let z[',num2str(i),'] := ',num2str(avgz(s)),'; \n'];
         end
         fprintf(fid,textv);    
         textv = [];
end %end avginitials


textv = [textv,' \n'];
textv = [textv,' \n'];
textv = [textv,' \n'];
textv = [textv,' \n'];
textv = [textv,' \n'];
textv = [textv,' \n'];
textv = [textv,' \n'];
textv = [textv,' \n'];

fprintf(fid,textv);
fclose(fid);

% tokens{1} = '#{n}#';
% tokens{2} = '#{nplus}#'; %there are n+1 ci's
% 
% replaces{1} = num2str(n);
% replaces{2} = num2str(n+1);
% 
% tokens{3} = '#{closeparam}#';
% if cc > 0
%     texty = ['#close contacts parameter. sum of close contacts constraint for target structure\n'];
%     texty = [texty,'param cc;\n'];
%     replaces{3} = sprintf([texty,'\n']);
% else
%     replaces{3} = '';
% end
% 
% tokens{4} = '#{bondls}#';
% %bond length constraints
% if bondls == 1
%     textz = ['## bond length constraints!! (very tricky - have to be built up from the first and second bonds)\n'];
%     textz = [textz,'var s{k in 1..',num2str(n-2),'} = sqrt(1-c[k]^2);\n'];
%     
%     textz = [textz,'var bu1{k in 1..',num2str(n-3),'} = tprime[k]*s[k]*s[k+1];\n'];
%     textz = [textz,'var bu2{k in 1..',num2str(n-3),'} = c[k]*c[k+1]-ti[k]*s[k]*s[k+1];\n'];
%     textz = [textz,'var bu3{k in 1..',num2str(n-3),'} = c[k+1];\n'];
%     textz = [textz,'\n'];
%     textz = [textz,'var bond{k in 1..',num2str(n-1),',j in 1..3};\n'];
%     textz = [textz,'\n'];
%     %textz = [textz,'#define the bond angles recursively: \n'];
%     %textz = [textz,'for {k in 1..',num2str(n-3),'} {   \n'];
%     textz = [textz,'var a1{k in 1..',num2str(n-3),'} = bond[k,2]*bond[k+1,3]-bond[k,3]*bond[k+1,2];\n'];
%     textz = [textz,'var a2{k in 1..',num2str(n-3),'} = bond[k,3]*bond[k+1,1]-bond[k,1]*bond[k+1,3];\n'];
%     textz = [textz,'var a3{k in 1..',num2str(n-3),'} = bond[k,1]*bond[k+1,2]-bond[k,2]*bond[k+1,1];\n'];
%     textz = [textz,'var deta{k in 1..',num2str(n-3),'} = a1[k]*bond[k,2]*bond[k+1,3] + a2[k]*bond[k,3]*bond[k+1,1] + a3[k]*bond[k,1]*bond[k+1,2]\n'];
%     textz = [textz,'           - a3[k]*bond[k,2]*bond[k+1,1] - a2[k]*bond[k,1]*bond[k+1,3] - a1[k]*bond[k,3]*bond[k+1,2];  \n'];
%     textz = [textz,'   \n'];
%     textz = [textz,'var ainv11{k in 1..',num2str(n-3),'} = (bond[k,2]*bond[k+1,3] - bond[k,3]*bond[k+1,2])/deta[k];   \n'];
%     textz = [textz,'var ainv21{k in 1..',num2str(n-3),'} = (bond[k+1,1]*bond[k,3] - bond[k,1]*bond[k+1,3])/deta[k];   \n'];
%     textz = [textz,'var ainv31{k in 1..',num2str(n-3),'} = (bond[k,1]*bond[k+1,2] - bond[k,2]*bond[k+1,1])/deta[k];   \n'];
%     textz = [textz,'   \n'];
%     textz = [textz,'var ainv12{k in 1..',num2str(n-3),'} = (a3[k]*bond[k+1,2] - a2[k]*bond[k+1,3])/deta[k];  \n'];
%     textz = [textz,'var ainv22{k in 1..',num2str(n-3),'} = (a1[k]*bond[k+1,3] - a3[k]*bond[k+1,1])/deta[k];  \n'];
%     textz = [textz,'var ainv32{k in 1..',num2str(n-3),'} = (a2[k]*bond[k+1,1] - a1[k]*bond[k+1,2])/deta[k];  \n'];
%     textz = [textz,'   \n'];
%     textz = [textz,'var ainv13{k in 1..',num2str(n-3),'} = (a2[k]*bond[k,3] - a3[k]*bond[k,2])/deta[k];  \n'];
%     textz = [textz,'var ainv23{k in 1..',num2str(n-3),'} = (a3[k]*bond[k,1] - a1[k]*bond[k,3])/deta[k];  \n'];
%     textz = [textz,'var ainv33{k in 1..',num2str(n-3),'} = (a1[k]*bond[k,2] - a2[k]*bond[k,1])/deta[k];  \n'];
%     textz = [textz,'   \n'];
%     textz = [textz,'var bn1{k in 1..',num2str(n-3),'} = ainv11[k]*bu1[k] + ainv12[k]*bu2[k] + ainv13[k]*bu3[k];   \n'];
%     textz = [textz,'var bn2{k in 1..',num2str(n-3),'} = ainv21[k]*bu1[k] + ainv22[k]*bu2[k] + ainv23[k]*bu3[k];   \n'];
%     textz = [textz,'var bn3{k in 1..',num2str(n-3),'} = ainv31[k]*bu1[k] + ainv32[k]*bu2[k] + ainv33[k]*bu3[k];   \n'];
%     textz = [textz,'   \n'];
%     textz = [textz,'#setting bnorm to 0 makes the program work, but I think it takes away the protection against\n'];
%     textz = [textz,'#error propogation. but better to be a bit wrong than to be nothing at all...\n'];
%     textz = [textz,'#the second and third columns of bond are a bit off without this, though \n'];
%     textz = [textz,'#var bnorm{k in 1..',num2str(n-3),'} >=.0001;# sqrt(bn1[k]*bn1[k] + bn2[k]*bn2[k] + bn3[k]*bn3[k]);     \n'];
%     textz = [textz,'     \n'];
%  
%     
%  %   textz = [textz,' } \n'];
%     textz = [textz,'\n'];
%     textz = [textz,'   \n'];
%     
%     replaces{4} = sprintf(textz);
% else %no bond constraints
%     replaces{4} ='';
% end
% 
% tokens{5} = '#{calphas}#';
% if cc > 0
%     texti = ['var calphas{i in 1..',num2str(n),', j in 1..3};\n'];  
%     texti = [texti,'#only sum up half, bc symmetric\n'];
%     texti = [texti,'var firstccsum{i in 1..',num2str(n),'} = sum{k in 1..i-1} 1/(.000001 + max(0.0001,(calphas[i,1]-calphas[k,1])^2 + (calphas[i,2]-calphas[k,2])^2 + (calphas[i,3]-calphas[k,3])^2));\n'];
%     texti = [texti,'#var otherccsum{i in 1..',num2str(n),'} = sum{k in i+1..',num2str(n),'} 1/(.000001 + max(.0001,(calphas[i,1]-calphas[k,1])^2 + (calphas[i,2]-calphas[k,2])^2 + (calphas[i,3]-calphas[k,3])^2));\n'];
%     texti = [texti,'var fccsum = sum{i in 1..',num2str(n),'} firstccsum[i];\n'];
%     texti = [texti,'#var occsum = sum{i in 1..',num2str(n),'} otherccsum[i];\n'];
%     
%     %texti = [texti, var ccconstraint = 
%     replaces{5} = sprintf(texti);
% else
%     replaces{5} = '';
% end
% 
% tokens{6} = '#{bondlconstraints}#';
% if bondls == 1
%     textz = [];
%     textz = [textz,'##bond length constraints\n'];
%     textz = [textz,'subject to bond11 : bond[1,1] == 1;\n'];
%     textz = [textz,'subject to bond12 : bond[1,2] == 0;\n'];
%     textz = [textz,'subject to bond13 : bond[1,3] == 0;\n'];
%     textz = [textz,'subject to bond21 : bond[2,1] == c[1];\n'];
%     textz = [textz,'subject to bond22 : bond[2,2] == s[1];\n'];
%     textz = [textz,'subject to bond23 : bond[2,3] == 0;\n'];
%     textz = [textz,'    \n'];
%     textz = [textz,'#defines the bond lengths   \n'];
%     textz = [textz,'subject to bonddef1{ k in 1..',num2str(n-3),'}:   bond[k+2,1] == bn1[k];  #/bnorm[k];    \n'];
%     textz = [textz,'subject to bonddef2{ k in 1..',num2str(n-3),'}:   bond[k+2,2] == bn2[k];  #/bnorm[k];      \n'];
%     textz = [textz,'subject to bonddef3{ k in 1..',num2str(n-3),'}:   bond[k+2,3] == bn3[k];  #/bnorm[k];      \n'];
%     textz = [textz,'    \n'];
%     
%     textz = [textz,'#restrictions on bond lengths: (so I can actually get rid of definitions of bonds)    \n'];
%     textz = [textz,'subject to lowerbondbound{ i in 1..',num2str(n-3),'}: bond[i + 2,1]^2 + bond[i + 2,2]^2 + bond[i + 2,3]^2  >= .95^2;     \n'];
%     textz = [textz,'subject to upperbondbound{ i in 1..',num2str(n-3),'}: bond[i + 2,1]^2 + bond[i + 2,2]^2 + bond[i + 2,3]^2  <= 1.05^2;     \n'];
%     textz = [textz,'# I think norm(bond(i+2)) is always 1, by the way geometry2bond was defined.     \n'];
%     textz = [textz,'    \n'];
%     
%     replaces{6} = sprintf(textz);
% else %no bond constraints
%     replaces{6} ='';
% end %tokens{16}
% 
% tokens{7} = '#{calphaconstraints}#';
% if cc > 0
%     textz = ['subject to calpa1{j in 1..3}: calphas[1,j] == 0;\n'];
%     textz = [textz,'subject to cfinder{k in 1..',num2str(n-1),',j in 1..3}: calphas[k+1,j] == calphas[k,j] + bond[k,j];\n'];
%     textz = [textz,'#close contact constraints:\n'];
%     textz = [textz,'#maybe it would be smarter to only sum of the first half of cc, too-- would\n'];
%     textz = [textz,'#probably save time in the long run -- okay, that is what i am doing now - thanks for the tip, me! \n'];
%     textz = [textz,'subject to ccupper: fccsum  <= cc + ',num2str(cc),';\n']; %   sum{i in 1..',num2str(n),'} \n'];
%     %textz = [textz,'           sum{k in 1..',num2str(n),'} ((calphas[i,1]-calphas[k,1])^2 + (calphas[i,2]-calphas[k,2])^2 + (calphas[i,3]-calphas[k,3])^2)^(-1/2) <= cc + 500;\n'];
%  
%     replaces{7} = sprintf(textz);
% else
%     replaces{7} = '';
% end
% 
% tokens{8} = '#{bondinitials}#';
% if bondls == 1
%     textz = ''; %normalize bond first -- bc geometry determines normed bonds
%     bond = bond./repmat(sqrt(sum(bond.^2,2)),1,3);
%     for i = 1:size(bond,1)
%         textz = [textz,'let bond[',num2str(i),',1] := ',num2str(bond(i,1)),';\n'];
%         textz = [textz,'let bond[',num2str(i),',2] := ',num2str(bond(i,2)),';\n'];
%         textz = [textz,'let bond[',num2str(i),',3] := ',num2str(bond(i,3)),';\n'];
%     end 
%     replaces{8} = sprintf(textz);
% else
%     replaces{8} = '';
% end
% 
% tokens{9} = '#{ccinitials}#';
% if cc == 1
%     textz = ''; %normalize bond first -- bc geometry determines normed bonds
%     bond = double(data{p}.bond);
%     bond = bond./repmat(sqrt(sum(bond.^2,2)),1,3);
%     x = calphas(bond);
%     for i = 1:size(x,1)
%         textz = [textz,'let calphas[',num2str(i),',1] := ',num2str(x(i,1)),';\n'];
%         textz = [textz,'let calphas[',num2str(i),',2] := ',num2str(x(i,2)),';\n'];
%         textz = [textz,'let calphas[',num2str(i),',3] := ',num2str(x(i,3)),';\n'];
%     end 
%     replaces{9} = sprintf(textz);
% else
%     replaces{9} = '';
% end
% 
% 
% tokens{10} = '#{initialc}#';
% tokens{11} = '#{initialt}#';

% else
%     replaces{10} = '';
%     replaces{11} = '';
% end
% 
% tokens{12} = '#{avginitialc}#';
% tokens{13} = '#{avginitialt}#';
% 
% if (avginitials)
%     %use average values of c, t, tp for each geo. class
%     load classavg;
%     avgc = classavgs(:,1);
%     avgcp = classavgs(:,2);
%     avgt = classavgs(:,3);
%     %find class sequence for this aa (already given as input for ampl in protein.dat, so not
%     %cheating)
%     claseq = classseq(data{p}.seq,lookup(:,:,:,:,end)); %class sequence for this p
%     sizeseq = size(claseq,1); %length of class sequence
%     
%     textv = ['# initial values. using avg values of c,cp,t,trpime for each class..\n'];
%     %%uh-oh! not so sure about averaging c and cplus 
%     i = 1;
%     g = claseq(i);
%     textv = [textv,'let c[',num2str(i),']:=',num2str(avgc(g),formatstring),';\n'];
%     for i = 2:sizeseq
%         g = claseq(i);
%         gb = claseq(i-1);
%         textv = [textv,'let c[',num2str(i),']:=',num2str((avgc(g)+avgcp(gb))/2,formatstring),';\n'];
%     end
%     i = sizeseq+1;
%     g = claseq(i-1);
%     textv = [textv,'let c[',num2str(i),']:=',num2str(avgcp(g),formatstring),';\n'];
%     replaces{11} = sprintf([textv,'\n']); % sprintf damit \n in echt umbrüche gewandelt werden
% 
%     %% set the ts to start very close to their original values
%     textv = [];
%     for i = 1:sizeseq
%         g = claseq(i);
%         textv = [textv,'let ti[',num2str(i),']:=',num2str(avgt(g),formatstring),';\n'];
%     end
%     replaces{12} = sprintf([textv,'\n']);
% 
%     %% and set the tprimes to start very close to their original values
%     textv = [];
%     for i = 1:sizeseq
%         g = claseq(i);
%         textv = [textv,'let tprime[',num2str(i),']:=',num2str(avgtprime(g),formatstring),';\n'];
%     end
%     replaces{13} = sprintf(textv);
% else
%     replaces{11} = '';
%     replaces{12} = '';
%     replaces{13} = '';
% end
% 
% 
% 
% 
% % write tokens
% global ampldir;
% 
% modfile = replaceinfile('modfile.mod.tpl',tokens,replaces,ampldir(1:end-1)); % remove trailing \
% 
% [pathstr, name, ext, versn] = fileparts(modfile);
% 
% tofile = [pathstr,filesep,problemname,'.mod'];
% 
% movefile(modfile,tofile);
% 
% dprintf('Renamed file %s.',tofile);
% 
% open(tofile);
% 
% 
