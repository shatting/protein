function [ bond ] = geometry2bond_ad( c, t, tp, len)
%GEOMETRY2BOND calculate normalized bond vectors from geometry c,t,t'
        %   first bond vector is parallel to x-axis, second lies in
        %   xy-plane.
        naa = length(c)+2;
        nfrags = naa-3;                                                     
        
        x = gradientinit([c;t;tp]);   
        % c(i) == x(i) (i \in 1..nfrags+1)
        % t(i) == x(nfrags+1+i) (i \in 1..nfrags)
        % tp(i) == x(2*nfrags+2+i) (i \in 1..nfrags)
        
        s = sqrt(1-x(1:naa-2).^2);
        
        A = gradient(zeros(3));
        bond = gradient(zeros(naa-1,3));
        bond(1,1) = 1;
        bond(2,1) = x(1);
        bond(2,2) = s(1); 
        
        for i=1:nfrags,            
            xci = i;
            xti = nfrags+1+i; % index of t_i in x
            xtpi = 2*nfrags+1+i;
            if (x(xti).x ~= t(i) || x(xtpi).x ~= tp(i) || x(xci).x ~= c(i))
                error('error');
            end
            
            A(1,:) = bond(i,:); % r_i
            A(2,:) = bond(i+1,:); % r_{i+1}
            A(3,:) = [bond(i,2)*bond(i+1,3)-bond(i,3)*bond(i+1,2)
                      bond(i,3)*bond(i+1,1)-bond(i,1)*bond(i+1,3)
                      bond(i,1)*bond(i+1,2)-bond(i,2)*bond(i+1,1)];            
            
            ss = s(i)*s(i+1);
%             b(1,1) = c(i)*c(i+1) - t(i)*ss;
%             b(2,1) = c(i+1);
%             b(3,1) = tp(i)*ss; %
            b(1,1) = x(xci)*x(xci+1) - x(xti)*ss;
            b(2,1) = x(xci+1);
            b(3,1) = x(xtpi)*ss; %
            
            bi = (ad_3x3inverse(A)*b)';   
            nbi = sqrt(sum(bi.*bi));            
            bond(i+2,:) = bi/nbi; % just to stabilize
            
        end
        
        for i=1:size(bond,1)
            bond(i,:) = bond(i,:)*len(i);
        end
end

