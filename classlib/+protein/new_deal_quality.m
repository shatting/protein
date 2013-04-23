function [ q, m ] = new_deal_quality( target, pred, print, maxpercent, sort )
%NEW_DEAL_QUALITY 

if (nargin<4), maxpercent = 100; end
if (nargin<5), sort = 0; end

if (sort)
    target = sortcl(target);
    pred = sortcl(pred);
end

nt = max(target);
np = max(pred);

conf = suffclass.utils.getfreq([target pred]);

if (nt == 4)
    disp('names:');
    disp(protein.utils.prednames(conf));
end

conf_max = conf./repmat(max(conf,[],1),max(target),1);

m = round(conf_max*100);

q = sum(conf_max ~= 1,2)/size(conf_max,2);
q2 = mean(conf_max,2);

if (nargin>2 && print)
    disp('frequencies:');
    suffclass.display.tightmat(freqs(target)');
    suffclass.display.tightmat(freqs(pred)');
    disp('confusion matrix:');
    suffclass.display.tightmat(conf);
    disp('conf_max (columns divided by max column entry):');
    suffclass.display.tightmat(m);
    dprintf('conf_max>=%i%%:',maxpercent);    
    starmat(m>=maxpercent);       
    
end

q = protein.kldiv.kldivqualitymeasure(target, pred, 0);

end
