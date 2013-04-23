function c = calphas(bondvector)
% convert from bond to calpha coordinates

bondvector = double(bondvector);

k=size(bondvector,1)+1;
c=zeros(k,3);

for i=2:k
    c(i,:)=c(i-1,:)+bondvector(i-1,:);
end
