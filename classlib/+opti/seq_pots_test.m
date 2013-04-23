function [ passed,v,vseq ] = seq_pots_test( optidata, detailed, version )
%SEQ_POTS_TEST Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    detailed = 0;
end

dprintf('- seq_pots_test.m');
Vseq_pots = {};
Vsuff_data2v = {};
totaltime = 0;
totalnfrags = data.GeomDB.db.nfrag;
for i=1:data.GeomDB.nch,
    chaini = data.GeomDB.db.chains{i}.getdatasetclone;
    [havefiltered, filteridx] = ...
        chaini.filter({'aa1','aa2','aa3','aa4'},@(feats)max(feats,[],2) <= 20);    
    if (havefiltered)
        totalnfrags = totalnfrags - sum(filteridx);
        if (detailed), dprintf('filtered %2i of %4i fragments in chain %4i',sum(filteridx),length(filteridx),i);end
        if ~chaini.getndata,
            if detailed, dprintf('..and skipped it (no aas <=20)'); end
            passed(i) = 1;
            continue;
        end
    end
    
    sgseq = optidata.sgclass.classifier.classify(chaini);
    sgcovseq = optidata.sgcovmodel(sgseq);
    
    chaini.ztransform(optidata.sclass.classifier,1);
    x = chaini.getknitrofeatures();
    
    switch (version)
        case 'fevaldata2'
            fevaldata = opti.FEvalData(sgcovseq);
            tic
            [v(i), Vseq_pots{i}] = opti.seq_pots_fevaldata(fevaldata,x);
            totaltime = totaltime + toc;  
        case 'fevaldata'
            fevaldata = opti.FEvalData(sgcovseq);
            tic
            [v(i), Vseq_pots{i}] = fevaldata.Vseq(x,0);
            totaltime = totaltime + toc;            
        case 'seq_pots_fast'            
            [ nfrags, idx, idx2, means, Ls, individualents ] = opti.seq_pots_fast_prepare( sgcovseq );
            tic
            [v(i), Vseq_pots{i}] = opti.seq_pots_fast(x, nfrags, idx, idx2, means, Ls, individualents );        
            totaltime = totaltime + toc;            
    end
    % compare to potential values from suffclass
    %[ cl, vmin, V, prob, Vdiff ] = suff_data2v( data, pot, pr, options )
    % use uniform prior and assume 1 potential per fragment
    data = chaini.getdata({'c','cp','z'});
    Vsuff_data2v{i} = zeros(1,length(sgseq));
    for j=1:length(sgseq)    
        pot = struct;
        pot.numpot = sgcovseq(j);
        pot.ncl = 1;
        pot.num = [1 2 3];
        pot.cat = [];    
        pot.typ = [0 0 0];
        [ cl, vmin, V] ...
            = suff_data2v( data(j,:), pot );        
        Vsuff_data2v{i}(j) = V;
    end
    %Vsuff_data2v{i} = diag(V);
        
    err = Vseq_pots{i} - Vsuff_data2v{i};
    passed(i) = max(abs(err)) < 1e-10;
    
    if (~passed(i) && detailed)
        str = sprintf('chain %i did not pass, maxabserr=%.2f',i,max(abs(err)));
        dprintf(str);
        figure;
        relerr = abs(err)./Vsuff_data2v{i};
        qot = Vseq_pots{i}./Vsuff_data2v{i};
        hist(relerr,50);
        plot(Vseq_pots{i},'r');
        hold on;
        plot(Vsuff_data2v{i},'b');
        plot(relerr,'k');
        plot(qot,'g');
        hold off;
        legend({'vseq fragment pot values','suff_data2v pot values','rel error','quotient'});
        title(str); 
    end
end

passed = all(passed);

if (detailed)
    dprintf('"%s" took %s, %i fragments per second.',version,showtime(totaltime),round(totalnfrags/totaltime));
end

if passed,    
    dprintf('PASSED.');
else
    dprintf('PROBLEMS FOUND!');
end

end

