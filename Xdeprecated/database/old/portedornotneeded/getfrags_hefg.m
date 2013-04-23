function [ features, torsold, cosold, data ] = getfrags_hefg( data )
%GEOM2FRAG Summary of this function goes here
%   Detailed explanation goes here

nump = length(data);

nf = 0;  % number of frags
for i=1:length(data),
    nf = nf + length(data{i}.seq)-3;
end

nfeats = 8;

dprintf('collecting fragment data, estimated %i fragments',nf);
% get geometry data
features = zeros(nf,nfeats,'int8');
torsold = zeros(nf,2);
cosold = zeros(nf,2);
count = 0;


for p=1:nump,

    if (mod(p,100) == 0),
        disp(sprintf('chain %i',p));
    end
    
    chainstruct = data{p};
    seq = chainstruct.seq;
    numaa = length(seq);
    
    %if ~isfield(chainstruct,'cosold'),
        chainstruct = new_deal_hefg_dataprep(chainstruct,0);
    %end

    data{p} = chainstruct;
    
    % get out of struct for speed
    cosa = chainstruct.cos;
    torsa = chainstruct.tors;
    torsolda = chainstruct.torsold;       
    cosolda = chainstruct.cosold;
    
    % get fragment primary sequences
    fragseq = zeros(numaa-3,4,'uint8');
    for i=1:numaa-3,
        fragseq(i,:) = seq(i:i+3);
    end
    
    maxfragaa = max(fragseq,[],2);
    idx = maxfragaa<=20;
    nidx = find(idx);
    
    if (nidx)
        nchfrag = double(sum(idx));
        % transpose into frag index
        fidx = count+1:count+nchfrag;
        
        features(fidx,1) = cosa(nidx);
        features(fidx,2) = cosa(nidx+1);
        features(fidx,3:4) = torsa(nidx,:);
        features(fidx,5:8) = fragseq(nidx,:);
        torsold(fidx,:) = torsolda(nidx,:);
        
        c1 = cosolda(1:end-1);
        c2 = cosolda(2:end);
        cosold(fidx,:) = [c1(nidx) c2(nidx)];        
        
        count = max(fidx);
    else
       dprintf('skipped chain %i, no fragments <= 20', p);
    end

end %end for p = 1:datasize

features(count+1:end,:) = [];
torsold(count+1:end,:) = [];
cosold(count+1:end,:) = [];