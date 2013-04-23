
if 0,
  if (exist('db')~=1),
    cd ..;startup;cd allesneu
    db = data2db( 'RT071127.mat' );    
  end
  if ~exist('fragdb'),
      fragdb = db2fragdb(db);
  end
elseif ~exist('fragdb'),
  tic
  disp('load data')
  load fragdb
  showtime
  %  nfrag: 221491
  %  freq(20,20,20,20)       frequency of primary sequence
  %  frags{20,20,20,20} 
  %      count               frequency
  %      data(count,nfeat)   features   
  %  nfeat: 6
  %  feats: {'cos(i)'  'cos(i+1)'  'tors(i,3)'  'len(i)'  'len(i+1)'  'res'}
  disp(fragdb)
end;

if ~exist('geom'),
    [ geom,seq ] = fragdb2linear( fragdb );
end
[HEL27] = HEL27_anal08(geom);

if 1,
    seq = double(seq);    
    featfun = @feat_pairs;
else
    %fseq = [seq(:,2:3), seq(:,1) + 20*(seq(:,4)-1)];
    featfun = @(x) x;
end
        
fseq = featfun(seq);

options = struct('plots',0,'max_iterations',20,'ask',0);
class_info = suff_classify(max(fseq,[],1),fseq,HEL27, options);
class_info.seq = seq; % save original fragments belonging to classes
class_info.featfun = featfun;

class_info.readme =  [
'        pot: [1x1 struct]      '
'         cl: [221491x6 double] '
'     confus: [27x27x6 double]  '
'       freq: [6x27 double]     '
'    options: [1x1 struct]      '
'        seq: [221491x4 uint8]  '
'    featfun: frag->features    '
'                               '];

save testclassify class_info

c = class_info.confus(:,:,end);
cp = c./repmat(max(c,[],1),27,1);
cs = c./repmat(max(c,[],2),1,27);
spy(cp>0.4);
title 'conf./repmat(max(conf,[],1),27,1) > 0.4';
figure;
spy(cs>0.8);
title 'conf./repmat(max(conf,[],2),1,27) > 0.8';



return;

data = seq;
target = HEL27;

[lookup,pot,confusion]=catclass([data HEL27], target);

if (~exist('lookupcluster'))
    dft = mkdft([data HEL27]);
    [lookupcluster,pairfreq,potcluster]=catcluster(dft,max(HEL27));
end

for i=1:size(data,1),
    
    v = [data(i,:), HEL27(i)];
    target(i) = lookupcluster.g(v(1),v(2),v(3),v(4),v(5));
    
end

[lookup,pot,confusion]=catclass([data HEL27], target);

