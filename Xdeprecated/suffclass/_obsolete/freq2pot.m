
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% freq2pot.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pot=freq2pot(pairfreq);
% create potential from a pair frequency table
%
% pairfreq(yi,yk,g,ik)  frequency of y with y_i=yi,y_k=yk in group g
%
% pot(yi,yk,g,ik)       potential for y with y_i=yi,y_k=yk in group g 
%                       ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
%                       i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
%                       k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
%
function pot=freq2pot(pairfreq);

ng=size(pairfreq,3);
on=ones(1,ng);
pairftot=sum(pairfreq,3);pairftot=pairftot(:,:,on,:);
pot=inf+zeros(size(pairfreq));
ind=(pairftot>0);
pot(ind)=-log(max(pairfreq(ind),realmin)./pairftot(ind));

