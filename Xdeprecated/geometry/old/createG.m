function [G] = createG(R)
n=size(R,1)+1;
G=zeros(n-1,4);

for i=1:n-4
   for j=1:4
     G(i,j)=R(i,:)*R(i+j-1,:)'; %%I want G to be an (n-1)x4 matrix
   end
end
  
for i=n-3:n-1
  G(i,1)=R(i,:)*R(i,:)';
end

G(n-3,2)=R(n-3,:)*R(n-2,:)';
G(n-3,3)=R(n-3,:)*R(n-1,:)';
G(n-2,2)=R(n-2,:)*R(n-1,:)';

end
