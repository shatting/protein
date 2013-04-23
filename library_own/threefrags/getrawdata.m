% get raw data information and save it for later use
% actual output depens on gramidx (which wont change i think), ninter and rbins(1)
% and rbins(end), since a fragment pair wont be looked at if the distance of
% the center aas isnt in latter interval.

ninter = 1; % minimum number of aa in between

%      1     6    11    16    21
%      2     7    12    17    22
%      3     8    13    18    23
%      4     9    14    19    24
%      5    10    15    20    25
gramidx = [6,11,12,16,17,18,21,22,23,24];

%rbins = 2:0.5:6 %~0.55% of all pairs, 170000
%rbins = 2:0.5:7 %~1.1% of all pairs, 340000, peak at 6.3
rbins = 2:0.5:10; %~4.08% of all pairs, 31190638, peaks at 6.3, low at 7, peak at 8.6, low at 9.2, peak at 10
%rbins = 5:0.5:10; %~3.8% of all pairs, 1,21mil, peaks at 6.3, low at 7, peak at 8.6, low at 9.2
nrbins = length(rbins)-1

nch = length(data);

npairs = 0;
ndatapairs = 0;

% suffs = cell(20,20,length(rbins)-1);
% emptysuff =  suffstat(1,zeros(10,1));

rawdata = struct;

nestimate = 148188;
dstore = zeros(nestimate,10+1); % gram(gramidx), rnorm
istore = uint8(zeros(nestimate, 7));    % seq(i,i+1,i+2,j,j+1,j+2), rnormbin
rawdata.readme = ['only threefrag pairs with rbins(1)<dist(i+1,j+1)<=rbins(end)     '
                  'dstore(l,:) = [gram(gramidx) dist(i+1,j+1)]  double              '
                  '              [   1:10             11     ]                      ' 
                  'istore(l,:) = [seq(i:i+2) seq(j:j+2) rbin ]  uint8               '
                  '              [    1:3        4:6      7  ]                      '
                  'generated with getrawdata.m                                      '];
    

for ch = 1:nch,
    if ~mod(ch,100), dprintf('chain %i, got %i datapoints (%2.2f%%)',ch, ndatapairs,100-(npairs-ndatapairs)/npairs*100); end
    bond = double(data{ch}.bond)/10;
    coords = bond2coords(bond);    
    
    seq = data{ch}.seq;
    naa = length(seq);
    %gram = bond*bond';
    
    
    for i=1:naa-5-ninter,
        
        for j=i+3+ninter:naa-2,
            npairs = npairs+1;
            r = coords(j+1,:) - coords(i+1,:);
            nr = norm(r);
            
            if nr<=rbins(1),
                %dprintf('dist(%i,%i)=%2.2f<=%2.1f in chain %i',i+1,j+1,nr,rbins(1),ch);
                continue;
            elseif nr>rbins(end) || seq(i+1)>20 || seq(j+1)>20,
                continue % not interested in these
            end
            
            ndatapairs = ndatapairs+1;
            % so r \in (rbins(1),rbins(end)]
                
            % get bin for r
            rbin = uint8(min(find(nr<=rbins))-1);
            
            % 
            % aa(i,i+1,i+2) ~ aa(j,j+1,j+2)
            % -> bond(i,i+1) ~ bond(j,j+1)
            % centers aa(i+1,j+1)
            % bonds i, i+1, j, j+1 and (coodrs(j+1)-coords(i+1))
            v = [bond(i:i+1,:); r; bond(j:j+1,:)];
            gram = v*v';
            feats = [gram(gramidx) nr];
                        
            if (ndatapairs>nestimate),
                disp('enlarging store');
                nestimate = ceil(nestimate*1.25);
                dstore(nestimate,1) = 0;
                istore(nestimate,1) = 0;
            end
            
            dstore(ndatapairs,:) = feats;
            istore(ndatapairs,:) = [seq([i:i+2,j:j+2])' rbin];
%             suff = suffs{seq(i+1),seq(j+1),rbin};
%             if (isempty(suff)),
%                 suff = emptysuff;
%             end
%             suff = suffstat(suff,feats,1,1);
%             
%             suffs{seq(i+1),seq(j+1),rbin} = suff;
        end        
        
    end
    

end

% truncate
if ndatapairs<nestimate,
   dstore(ndatapairs+1:end,:) = [];
   istore(ndatapairs+1:end,:) = [];
end

rawdata.dstore = dstore;
rawdata.istore = istore;
rawdata.rbins = rbins;
rawdata.gramidx = gramidx;
rawdata.npairs = npairs;

disp('saving');
save rawdata

rawdata