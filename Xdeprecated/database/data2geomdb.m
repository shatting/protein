function [ db ] = data2geomdb( rtfile )
% [ db ] = data2geomdb( rtfile )
% augment data with geometry data
% uses: bond2geomstruct.m

tic
load(rtfile);
fieldnames(data{1})
nch=length(data)
   
% simple statistics
naa=zeros(1,nch);
res=zeros(1,nch);
svdbb=zeros(3,nch);
afreq=zeros(24,1);
res24=inf;
for ch=1:nch,
  if rem(ch,100)==0,
    s = showtime;
    dprintf('computing geometry for chain %i (elapsed %s)',ch,s);
  end;
  dat=data{ch};

  % amino acid frequencies
  aa=dat.seq;% amino acid sequence
  naa(ch)=length(aa);
  res(ch)=dat.res; 
  % if res(ch)==0, res(ch)=inf;data{ch}.res=inf; end;
  for t=1:naa(ch)-3,
    afreq(aa(t))=afreq(aa(t))+1;    
  end;
  if max(aa)>23, res24=min(res24,res(ch)); end;

  data{ch} = bond2geomstruct(dat);
  
  svdbb(:,ch) = data{ch}.svdbb;
  
end;

readme=[
'% generated with data2geomdb.m                                 '
'% .geomdata{ch} contains (some of) the fields                  '
'%   .name    % name of subchain                                '
'%   .res     % resolution (-1 = NMR)                           '
'%   .seq     % amino acid sequence, uint8                      '
'%   .bond    % bond vectors, 3 columns, int8                   '
'%   and all fields from bond2geomstruct.m                      '
];

db.geomstructs = data;
db.nch = nch;
db.afreq = afreq;
db.readme = readme;
db.naa = sum(naa);
db.nfrag = sum(naa-3);
db.res24 = res24;
db.svdbb = svdbb;
db.geomversion = bond2geomstruct();

end

