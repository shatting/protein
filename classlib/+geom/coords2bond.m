function b = coords2bond(coords)

    k=size(coords,1)-1;
    b=zeros(k,3);

    for i=1:k
        b(i,:)=coords(i+1,:) - coords(i,:);
    end

end

