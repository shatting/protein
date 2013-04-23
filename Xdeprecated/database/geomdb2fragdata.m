function [ fragdata ] = geomdb2fragdata( geomdb, featurenames )
% [ fragdata ] = geomdb2features( geomdb, featfun )
% extract fragment data from geomdb.
% INPUT:    geomdb
%           featurenames       cell array of feature names or a
%                              featurefunction
%                           

dprintf('collecting fragment features from %i fragments',geomdb.nfrag);
% get geometry data

if (isempty(featurenames) || nargin < 2)
    % get field names
    [x,featurenames] = geomstruct_getfeaturesbyname(geomdb.geomdata{1});
end

if ~iscell(featurenames),
    featfun = featurenames;
else    
    featfun = @(geomstruct,idx)geomstruct_getfeaturesbyname(geomstruct,featurenames,idx);
end

features = zeros(geomdb.nfrag,length(featurenames)); 

seq = zeros(geomdb.nfrag,4);
count = 0;

for ch=1:geomdb.nch,

    if (mod(ch,100) == 0),
        disp(sprintf('collecting features from chain %i',ch));
    end
    
    geomstruct = geomdb.geomdata{ch};
    chseq = geomstruct.seq;    
        
    % get fragment primary sequences
    fragseq = seq2fragseq(chseq);
    
    maxfragaa = max(fragseq,[],2);
    idx = maxfragaa<=20;
    nidx = find(idx); % 
    
    if (nidx)
        nchfrag = double(sum(idx));
        % transpose into frag index
        fidx = count+1:count+nchfrag;
        
        features(fidx,:) = featfun(geomstruct,nidx);
        seq(fidx,:) = fragseq(nidx,:);
               
    else
       dprintf('skipped chain %i, no fragments <= 20', ch);
    end
    count = max(fidx);
    
end %end for ch = 1:datasize

% truncate
features(max(fidx)+1:end,:) = [];
seq(max(fidx)+1:end,:) = [];

fragdata = struct;
fragdata.data = features;
fragdata.seq = seq;
fragdata.featurenames = featurenames;

end

