

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getfreq.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function freq=getfreq(list);
% calculate frequency table from a list of examples
% freq(a,b,...) is the number of l with list(l,1)=a, list(l,2)=b, ...
%
% function freq=getfreq(list,wt);
% calculate frequency table from a weighted list of examples
% freq(a,b,...) is the sum of wt(i) with list(l,1)=a, list(l,2)=b, ...
%
function freq=getfreq(list,wt);


if min(list(:))<1, 
  minentry=min(list(:))
  error('entries must be positive integers'); 
end;

list=double(list);
[N,ndim]=size(list);
if nargin==1, wt=ones(N,1); end;

dims=max(list,[],1);
dim=prod(dims);
if dim>=2^31, 
  dims,product=dim
  error('number of cases must be <2^31'); 
end;

code=list(:,ndim);
for k=ndim-1:-1:1,
  code=(code-1)*dims(k)+list(:,k);
end;
freq=sum(sparse(1:N,code,wt,N,dim),1);
if (length(dims)>1)
    freq=reshape(full(freq),dims); % full needed if more than 2 indices
else
    freq=full(freq);
end


