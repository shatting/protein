function [ fragdb ] = db2fragdb( db )
%GEOM2FRAG Summary of this function goes here
%   Detailed explanation goes here

nump = db.nch;

disp(sprintf('collecting fragment data.'));
INITIAL_SIZE = 10;
nfeats = 6;

skip = 0;

fragstruct_proto.count = 0;
fragstruct_proto.size = INITIAL_SIZE;
fragstruct_proto.data = zeros(INITIAL_SIZE,nfeats);

% get geometry data

fourfrags = cell(20,20,20,20);
fourfrags_freq = zeros(20,20,20,20);
features = zeros(1,nfeats);

for p=1:nump,%(1:1334)

    if (mod(p,100) == 0),
        disp(sprintf('chain %i',p));
    end

    geom = db.geom{p};
    seq = geom.seq;
    numaa = length(seq);


    for i=1:numaa-3,
        aa1 = seq(i);
        aa2 = seq(i+1);
        aa3 = seq(i+2);
        aa4 = seq(i+3);

        %skip fragments with nonrealistic aa's
        if aa1 > 20 || aa2 > 20 || aa3 > 20 || aa4 > 20,
            skip = skip + 1;
            continue;
        end

        features(1)=geom.cos(i);
        features(2)=geom.cos(i+1);
        features(3)=geom.tors(i,3);
        features(4)=geom.len(i);
        features(5)=geom.len(i+1);
        features(6)=geom.res;
        
        fragstruct = fourfrags{aa1,aa2,aa3,aa4};
        
        %check if this combo of aa's has occurred yet:
        if ~isstruct(fragstruct), %if hasn't occurred:
            fragstruct = fragstruct_proto; %gives fragstruct structure properties: .count, .data, and .size
        else %if already there
            % make bigger?
            if fragstruct.count == fragstruct.size,
                fragstruct.size = fragstruct.size * 2;
                fragstruct.data(fragstruct.size,1) = 0;
            end
        end

        fragstruct.count = fragstruct.count + 1;
        fragstruct.data(fragstruct.count,:) = features;

        fourfrags{aa1,aa2,aa3,aa4} = fragstruct;

    end
end %end for p = 1:datasize

% truncate
for aa1=1:20,    

    for aa2=1:20,

        for aa3=1:20,
            for aa4=1:20,

                fragstruct = fourfrags{aa1,aa2,aa3,aa4};
                
                if isstruct(fragstruct),
                    
                    fourfrags_freq(aa1,aa2,aa3,aa4) = fragstruct.count;
                    
                    % truncate buckets
                    fragdata = fragstruct.data(1:fragstruct.count,:);
                    
                    fourfrags{aa1,aa2,aa3,aa4} = rmfield(fourfrags{aa1,aa2,aa3,aa4},'size');
                    fourfrags{aa1,aa2,aa3,aa4}.data = int16(fragdata);
                    fragstruct = [];
                end
                
            end
        end
    end
    disp(sprintf('%i %%, fillrate %f',round(100*aa1/20), sum(sum(sum(fourfrags_freq(aa1,:,:,:) ~= 0)))/20^3 ));
end

fragdb.nfrag = db.nfrag-skip;
fragdb.freq = fourfrags_freq;
fragdb.frags = fourfrags;
fragdb.nfeat = nfeats;
fragdb.feats = {'.cos(i)','.cos(i+1)','.tors(i,3)','.len(i)','.len(i+1)','.res'};
fragdb.readme = 'generated with db2fragdb.m. field names in .feats refer to fields generated with data2db.m';