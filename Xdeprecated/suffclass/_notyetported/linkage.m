
% written by Waltraud Huyer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% linkage.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ord,c,bounds,diss] = linkage(d,opt)
% hierarchical cluster analysis by single or complete linkage 
% 
% Input:
% d		symmetric matrix of dissimilarities between classes i 
%		and j
% opt 		opt=0: single linkage, opt=1: complete linkage
%
% Output:
% ord		permutation of the set of classes such that the fusions 
%		in the clustering algorithm can be performed without 
%		crossing lines
% c		In the k-th step of the cluster analysis, two clusters
% 		with dissimilarity c(k) are joined
% bounds	ord(bounds(1,k):bounds(2,k)) gives the elements of the
%		cluster generated in the k-th step
% diss		vector of length n-1 of dissimilarities between 
%		joining clusters; joining classes ord(j), ord(j+1) 
%		when diss(j) <= c(i) gives n - i clusters
%		c = sort(diss)
%		
function [ord,c,bounds,diss] = complink1(d,opt)
n = size(d,1);
for i=1:n
  d(i,i) = Inf;
end
ind0 = 1:n;
ord = 1:n;
diss = zeros(1,n-1);
len = ones(1,n);
pt = 1:n;
for k = 1:n-1
  [y,i] = min(d);
  [y,j] = min(y);
  i = i(j);
  if i < j
    i1 = j;
    j = i;
    i = i1;
  end
  if pt(i) > pt(j)+len(j) 
    [ord,diss,len,pt] = swap1(ord,diss,y,len,pt,j,i);
  else
    j1 = pt(j);
    diss(j1+len(j)-1) = y;
    len(j) = len(j) + len(i);
  end
  len1(k) = len(j);
  ind0(find(ind0==i)) = [];
  if length(ind0) > 1
    ind = ind0;
    ind(find(ind==j)) = [];
    if opt, d(j,ind) = max(d(j,ind),d(i,ind));
    else  d(j,ind) = min(d(j,ind),d(i,ind));
    end;
    d(ind,j) = d(j,ind)';
    d(i,1:n) = Inf*ones(1,n);
    d(1:n,i) = d(i,1:n)';
  end
end
[c,ind] = sort(diss);
for k = 1:n-1
  if k > 1
    j = find(bounds(2,:)==ind(k));
  end
  if k == 1 | isempty(j)
    bounds(:,k) = [ind(k); ind(k)+len1(k)-1];
  else
    j = j(length(j));
    bounds(:,k) = [bounds(1,j); bounds(1,j)+len1(k)-1];
  end
end





  
