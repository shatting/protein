% SUFF_DATA2CL Non-vectorized class prediction.
%
% clpred = SUFF_DATA2CL( data, pot, options )
%
% Calculate predicted class assignments for data under potential pot. Only
% works fast for categorical features with a lookup table.
% Otherwise, iterates suff_vectorpotential.
% For a vectorized version, see suff_data2v.
%
% INPUT:
%   data        nxd array, n feature vectors of dimension d
%   pot         potential struct (see suff_pot)
%   options     
%
% See also: suff_data2v, suff_pot, suff_vectorpotential.

function clpred = suff_data2cl( data, pot, options )

if (nargin < 3), options = []; end;    

[N,d]=size(data);

dprintf('computing new class assignments: %i items',N);
if (isempty(pot.num) && isfield(pot.catpot,'lookup')),
    datacat = data(:,pot.catpot);
    lookup = pot.catpot.lookup;

    [N,d]=size(datacat);
    if isa(datacat,'double'),
      clpred=lookup.g(datacat*lookup.fac-lookup.shift);
    else
      % convert to double in batches
      clpred=uint8(zeros(N,1));
      batchsize=100000;
      nn=fix(batchsize/d);
      nb=fix(N/nn)+1;
      for ind1=1:nn:N,
        ind=[ind1:min(N,ind1+nn-1)];
        clpred(ind)=lookup.g(double(datacat(ind,:))*lookup.fac-lookup.shift);
      end;
    end
else
    clpred = zeros(N,1);
    rep = fix(N/10);
    for i=1:N,
        if rem(i,rep)==0 && N>5000,
          %timetogo=showtime((toc-tim1)*togo/(l-1));
          dprintf('%.2f%%',i/N*100);
        end;
        [x, clpred(i)] = suff_vecpot(data(i,:), pot, options);
    end
end
% removed because of unlabelled items
%
%if nargout==1, return; end;
% % compute confusion matrix
% confusion=getconf(oldcl,clpred);
% if nargout==2, return; end;
% 
% % compute target prediction table
% if nargin<4,cost=1;end;
% if length(cost)==1,
%   [n,nn]=size(confusion);
%   if n~=nn, error('default size requires group=class'); end;
%   cost=1-eye(n);
% end;
% [mc,g2t]=min(cost*confusion);g2t=g2t(:)';
