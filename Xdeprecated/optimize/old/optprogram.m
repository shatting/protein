function [allpres,GDTs,bestsGDT,bestpreGDT,bestsRMSE,bestpreRMSE,S,sactual] = optprogram(p,confus,newclasses,class_info,data,probname)
%[GDTs,bestsGDT,bestpreGDT,bestsRMSE,bestpreRMSE,S,sactual] = optprogram(p,confus,newclasses,class_info,data,probname)
% finds likely class sequences for the protein p, 
% then optimizes with each of the possible sequences,
% and then chooses the sequences with the lowest rmsd and highest GDT
% as the winning sequence
% plots all structures
% returns rmse of all structures, the sequence (bests) with the best rmse,
% and the ti, ci values of the best sequence (bestpre.c, bestpre.t)

if newclasses == 1
    nums = 4;
else
    nums = 27;
end

numsequences = 19;  

% parameters for runopt
algnum = 1; maxit = 3000; clco = 1; hp = 1; initials = 0; avginitials = 1; thfrags = 0; zlu = 1; individualopt = 0; displaymode = 0;

%reset _boundparams.dat - just in case I was running boundchecker and
%stopped early and forgot
load sbounds;
global local_ampldir;
disp('restoring boundparams.dat to original --> the file for c,cp, and z boundaries');
individualopt = 0;
makeboundparamsdotdat([local_ampldir,filesep,'_boundparams','.dat'],nums,individualopt,clower,cupper,cplower,cpupper,zlower,zupper);

% need to store info from onesequenceforp, so that this doesn't waste
% calculations
[Smore,SVmore,CL,jsofeachs,j,i,pastis,CLcol,jeq1,L] = onesequenceforp(p,data,confus,class_info,newclasses,nums,numsequences);
% to get moresequences, run this later: 
%moresequences = 20
%[S,SV,CL,jsofeachs,j,i,pastis,CLcol,jeq1,L] = onesequenceforp(p,data,confus,class_info,newclasses,nums,moresequences,CL,jsofeachs,j,i,pastis,CLcol,Smore,SVmore,jeq1,L)
% doesn't work yet if original nfound is smaller than 3, or something

sv = SVmore;
S = Smore;


%old = double(data{p}.bond);
%old = old./repmat(sqrt(sum(old.^2,2)),1,3);
old = rotatematrix(data{p}.bond);

figure(1);
subplot(5,4,1); 
ribbon3d(calphas(old),5); 
title(plabel(0,sv(1))); %findV(p,data,confus,class_info,newclasses)));

%now plot prediction with actual sequence
geo = bond2geo(data{p}.bond);
sactual = geo.tcl;
runopt(data,class_info,newclasses,confus,probname,p,algnum,maxit,clco,hp,initials,avginitials,thfrags,zlu,displaymode,individualopt,sactual);  
   prez = fromampl(probname,sactual,newclasses); 
    if newclasses == 0    
        pregeom = geometryz2geometry(prez);   
        preactual = geometry2bondn(pregeom);
    else
        prez.tqu = prez.zo;
        prez.tcl = sactual;
        prez.co = prez.co(1:end-1);
        prez.cpo = prez.co(end);
        preactual = geo2bond(prez);
    end
    %[ Q, q, rmse ] = register(old,preactual)
    [GDT,LGArmse,lgacoords] = dolga(preactual,data{p});
figure(2);
subplot(1,2,1); 
ribbon3d(lgacoords,5); 
title('actual structure of protein vs predicted structure using actual geometry sequence (with lga superposition)'); %findV(p,data,confus,class_info)));
subplot(1,2,2);
ribbon3d(calphas(preactual),5);
title(plabel(LGArmse,0,GDT));

rmsemin = inf;
GDTmax = -inf;

bestsRMSE = 0; % to keep track of which sequence has best prediction 
           % according to rmse
bestpreRMSE = 0;
bestsGDT = 0; % to keep track of which sequence has best prediction 
           % according to GDT
bestpreGDT = 0;

rmses = zeros(1,numsequences); 
GDTs = zeros(1,numsequences); 

figure(4)
subplot(5,4,1)
title('lgas superposition of the template');


allpres = cell(numsequences,1);

for i = 1:numsequences
    i
    % input arguments at end of runopt: p,algnum,maxit,clco,hp,initials,zlu,s
    runopt(data,class_info,newclasses,confus,probname,p,algnum,maxit,clco,hp,initials,avginitials,thfrags,zlu,displaymode,individualopt,S(i,:));
    prez = fromampl(probname,S(i,:),newclasses); 
    if newclasses == 0    
        pregeom = geometryz2geometry(prez);   
        pre = geometry2bondn(pregeom);
    else
        prez.tqu = prez.zo;
        prez.tcl = S(i,:);
        prez.co = prez.co(1:end-1);
        prez.cpo = prez.co(2:end);
        pre = geo2bond(prez);
    end
    %geometry2bondn(fromampl(probname),S(i,:),newclasses);%*3.8;
%     disp('twist with Q and shift by q to "overlap" with original configuration:')
%     [ Q, q, rmse ] = register(old,pre)
   
    allpres{i} = pre;

    [GDT,LGArmse,lgacoords] = dolga(pre,data{p});
    
    rmse = LGArmse;
    
    figure(1)
    subplot(5,4,i+1), 
    ribbon3d(calphas(pre),5);  
    GDTs(i) = GDT;
    rmses(i) = LGArmse;
    title(plabel(rmse,sv(i),GDT));
    drawnow;
    figure(4)
    subplot(5,4,i+1), 
    ribbon3d(lgacoords,5); 
    title(plabel(rmse,sv(i),GDT));
    drawnow;
    
    if rmse < rmsemin %check if best rmse
        rmsemin = rmse;
        bestsRMSE = S(i,:);
        bestpreRMSE = pre;
        bestsvRMSE = sv(i);
        bestlgacoordsRMSE = lgacoords;
        bestpreRMSEGDT = GDT;
    end
    if GDT > GDTmax
        GDTmax = GDT;
        bestsGDT = S(i,:);
        bestpreGDT = pre;
        bestsvGDT = sv(i);
        bestlgacoordsGDT = lgacoords;
        bestpreGDTrmse = rmse;
    end
end
   
disp('best prediction (according to rmse) made by: ')
bestsRMSE
%makepdbform(1,3.8*calphas(bestpreRMSE),3.8*calphas(old),data{p}.seq);
figure(3)
subplot(1,3,1); 
ribbon3d(bestlgacoordsRMSE,5); 
title('best prediction by rmse, lga superposition of template onto model');
subplot(1,3,2);
ribbon3d(calphas(bestpreRMSE),5);
title(plabel(rmsemin,bestsvRMSE,bestpreRMSEGDT));

disp('best prediction (according to GDT) made by: ')
bestsGDT
%makepdbform(1,3.8*calphas(bestpreRMSE),3.8*calphas(old),data{p}.seq);
figure(5)
subplot(1,2,1); 
ribbon3d(bestlgacoordsGDT,5); 
title('best prediction by gdt, lga superposition of template onto model');
%title(plabel(0,bestsv)); %findV(p,data,confus,class_info)));
subplot(1,2,2);
ribbon3d(calphas(bestpreGDT),5);
title(plabel(bestpreGDTrmse,bestsvGDT,GDTmax));

