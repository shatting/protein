function [ allpres, mrmsds, actualsprob, HMMSeqs] = predict( pchain, optidata, probname, options)
% opti.predict( pchain, optidata, probname, options)
%   Rewrite of optprogram.m which did the following:predict( pchain, optidata, probname, options)
%
%[GDTs,bestsGDT,bestpreGDT,bestsRMSE,bestpreRMSE,HMMSeqs,sactual] = optprogram(proteinnumber,confus,newclasses,new_deal_data,database,problemname)
% finds likely class sequences for the protein proteinnumber, 
% then optimizes with each of the possible sequences,
% and then chooses the sequences with the lowest rmsd and highest GDT
% as the winning sequence
% plots all structures
% returns rmse of all structures, the sequence (bests) with the best rmse,
% % and the ti, ci values of the best sequence (bestpre.c, bestpre.t)
%
% 
% if ((proteinnumber == 98) || (proteinnumber == 99) || (proteinnumber == 100)) %these proteins have more than 20% aa's greater than 20
%     disp('this is a ridiculous protein, filled with unrealistic amino acids');    
%     return;
% end
% parameters for runopt
%algnum = 1; maxit = 3000; clco = 1; hp = 1; initials = 0; avginitials = 1; thfrags = 0; zlu = 1; individualopt = 0; displaymode = 0;
tic

if (options.fastobjective && options.savexhistory)
    dprintf('fastobjective and savexhistroy are not supported together. resorting to fastobjective only.');
end

naa = length(pchain.seq);
% need to store info from onesequenceforp, so that this doesn't waste
% calculations
chaingammaclass = class.Classification(pchain,optidata.gclass.classifier);
chainactualsclass = class.Classification(pchain,optidata.sclass.classifier);


HMMSeqs = [];
if options.HMM_numseqs > 0
    dprintf('%s: HMM: Calculating %i best s-class sequences',showtime,options.HMM_numseqs);
[HMMSeqs,HMMPots,CL,jsofeachs,j,i,pastis,CLcol,jeq1,L] = ...
    hmm.getnbest(chaingammaclass.cl,optidata.sgconfus,options.HMM_numseqs);    
end
% to get moresequences, run this later: 
%moresequences = 20
%[HMMSeqs,SV,CL,jsofeachs,j,i,pastis,CLcol,jeq1,L] = onesequenceforp(proteinnumber,database,confus,new_deal_data,newclasses,nums,moresequences,CL,jsofeachs,j,i,pastis,CLcol,HMMSeqs,HMMSeqs,jeq1,L)
% doesn't work yet if original nfound is smaller than 3, or something

% --------- RUN optimization with actual s class sequence

%chainsgclass = class.Classification(pchain,optidata.sgclass.classifier);
actualsprob = [];
if (options.sactual)
    dprintf('%s: Setting up actual s-class sequence knitro problem.',showtime);
    actualoptions = options;
    actualoptions.savexhistory = options.showsactualrmsdprogression;
    etime = 8*naa^3/1.9048e+004;
    %dprintf('%s: Running optimization on actual s-class sequence.',showtime);
    %dprintf('%s: Estimated time for one optimization run: %s',showtime,showtime(etime));
    actualsprob = opti.KnitroProblem(pchain.seq,chainactualsclass.cl,chaingammaclass.cl,optidata,actualoptions);
    actualsprob.run(options.getknitrotime);

    rrmsd = actualsprob.predcoords.registerto(pchain);
    %armsd = actualsprob.predcoords.rmsdto(pchain);
    dprintf('%s: min rmsd using actual s-class seq: %.2f Ang',showtime,rrmsd/10);
    %[GDT,LGArmse,lgacoords] = dolga(preactual,database{proteinnumber});
    figure(1);
    subplot(1+options.showsactualrmsdprogression,2,1); 
    pchain.draw;
    title('chain geometry');
    subplot(1+options.showsactualrmsdprogression,2,2); 
    actualsprob.predcoords.draw;
    title(sprintf('prediction (actual s classes).\nminimum rmse=%.2f',rrmsd/10));
    drawnow;
    if (options.showsactualrmsdprogression)
        subplot(2,2,[3 4]);    
        sactualrmsdprogression = actualsprob.getxhistoryrmsdsto(pchain);
        plot(sactualrmsdprogression);
        title('rmsd progression');
        drawnow;
    end
end
% ------------- RUN optimization for all HMM s sequences
% rmsemin = inf;
% GDTmax = -inf;
% 
% bestsRMSE = 0; % to keep track of which sequence has best prediction according to rmse
% bestpreRMSE = 0;
% bestsGDT = 0; % to keep track of which sequence has best prediction according to GDT
% bestpreGDT = 0;
% 
% rmses = zeros(1,options.HMM_numseqs); 
% GDTs = zeros(1,options.HMM_numseqs); 

%figure(4)
%subplot(5,4,1)
%title('lgas superposition of the template');

mrmsds = zeros(options.HMM_numseqs,1);
%armsds = zeros(options.HMM_numseqs,1);

dprintf('%s. preparing knitro problems.',showtime);
for i=1:options.HMM_numseqs,
    allpres(i) = opti.KnitroProblem(pchain.seq,HMMSeqs(i,:)',chaingammaclass.cl,optidata,options);
end
dprintf('%s. done.',showtime);

% determine parallel use
aminpool = matlabpool('size') > 0;
if (~aminpool)    
    dprintf('%s. Using single core.',showtime);    
else 
    dprintf('%s. Using %i cores.',showtime,matlabpool('size'));        
end

% main loop. catch exceptions to close matlabpool.
try
    tic;
    parfor i = 1:options.HMM_numseqs
        
        % set matlab process priority on each worker
        %if (~strcmp(options.priority,priority)), priority(options.priority); end
                
        if ~aminpool
            wname = '';
        elseif (isempty(getCurrentWorker))
            wname = 'worker';
        else
            wname = [get(getCurrentWorker,'Name'),': '];            
        end
        
        %dprintf('%s[HMM %3i] starting optimization.',wname,i);
        allpres(i) = allpres(i).run(options.getknitrotime);        
                
        mrmsds(i) = allpres(i).predcoords.registerto(pchain)/10; %Angstrom
        %test
        %xco = allpres{i}.predcoords.coords;    
        %armsds(i) = geom.Coords(geom.coords2bond(xco(end:-1:1,:))).registerto(pchain);
                
        dprintf('%s[HMM %3i] done. min rmsd=%3.2f\texitflag=%4i\tfirstorderopt=%g',wname,i,mrmsds(i),allpres(i).EXITFLAG,allpres(i).OUTPUT.firstorderopt);           
        
        
        %[GDT,LGArmse,lgacoords] = dolga(pre,geomdb.geomdata{proteinnumber});
                
        if (options.plotall)
            figure(i+1);
            subplot(1,2,1); 
            pchain.draw;
            title('chain geometry');
            subplot(1,2,2); 
            allpres(i).predcoords.draw;
            title(sprintf('prediction with HMM seq #%i.\nmin rmsd=%.2f Ang',i,mrmsds(i)));
        end
    end % parfor
    
    dprintf('%s: done.',showtime);    
       
catch e
    disp(e.getReport);
    dprintf('%s: error.',showtime);    
    if (matlabpool('size') > 0)
        matlabpool close;
    end
end

if options.HMM_numseqs>0
    [~,minidx] = min(mrmsds);
    str = sprintf('best prediction according to min rmsd (%.2f): HMM seq #%i.',mrmsds(minidx),minidx);
    
    if isfield(options,'plotresult') && options.plotresult
        figure;
        subplot(1,3,1); 
        pchain.draw;
        title('chain geometry');
        subplot(1,3,2);
        allpres(minidx).predcoords.draw;
        
        title(str);        
    end
    
    dprintf(str);
end
    
% disp('best prediction (according to GDT) made by: ')
% bestsGDT
% %makepdbform(1,3.8*calphas(bestpreRMSE),3.8*calphas(old),database{proteinnumber}.seq);
% figure(5)
% subplot(1,2,1); 
% ribbon3d(bestlgacoordsGDT,5); 
% title('best prediction by gdt, lga superposition of template onto model');
% %title(plabel(0,bestsv)); %findV(proteinnumber,database,confus,new_deal_data)));
% subplot(1,2,2);
% ribbon3d(calphas(bestpreGDT),5);
% title(plabel(bestpreGDTrmse,bestsvGDT,GDTmax));

end
