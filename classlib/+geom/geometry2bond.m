function [ bond, Acond ] = geometry2bond( c, t, tp, len, useintlabAD )
%GEOMETRY2BOND calculate normalized bond vectors from geometry c,t,t'
        %   first bond vector is parallel to x-axis, second lies in
        %   xy-plane.
% [ bond, Acond ] = geometry2bond( c, t, tp, len )        
% [ bond ] = geometry2bond( cgrad, tgrad, tpgrad, len )
    
    verifyhandAD = 0; % this is just useful for coding
    
    if (nargin < 5)
        useintlabAD = 0;
    end
    
    if (~isa(c,'gradient')) % [ bond, Acond ] = geometry2bond( c, t, tp, len )        
        naa = length(c)+2;
        nfrags = naa-3;
        
        s = sqrt(1-c.^2);        
        ss = s(1:end-1).*s(2:end);
        
        bond = zeros(naa-1,3);                
        Acond = zeros(naa-3,1);                
        
        bond(1,1) = 1;
        bond(2,1) = c(1);
        bond(2,2) = s(1); % TODO: *sign(c(1)); removed for consistency with deriving version
        
        for i=1:nfrags,            
            
            A(1,:) = bond(i,:); % r_i
            A(2,:) = bond(i+1,:); % r_{i+1}
            % A(3,:) = cross(bond(i,:),bond(i+1,:)); % r_i x r_{i+1}
            % above is awfully slow
            A(3,:) = [bond(i,2)*bond(i+1,3)-bond(i,3)*bond(i+1,2)
                      bond(i,3)*bond(i+1,1)-bond(i,1)*bond(i+1,3)
                      bond(i,1)*bond(i+1,2)-bond(i,2)*bond(i+1,1)];            
                        
            b(1,1) = c(i)*c(i+1) - t(i)*ss(i);
            b(2,1) = c(i+1);
            b(3,1) = tp(i)*ss(i);
            
            bond(i+2,:) = (A\b)';   
            bond(i+2,:) = bond(i+2,:)./norm(bond(i+2,:)); % stabilize
            
%             % this fixes a bad bug            
%             dets(i) = det(bond(i:i+2,:));
%             if (dets(i)<0 || 0)
%                 bond(i+2,[2,3]) = -bond(i+2,[2,3]);
%             end
            if (nargout > 1)
                Acond(i) = cond(A);
            end
        end
        % not necessary anymore after stabilization was introduced
%         lengthdeviations = abs(1-sqrt(sum(bond.^2,2)));
%         maxdeviation = max(lengthdeviations);
%         if maxdeviation > 1e-8
%             dprintf('reconstructed bond length has maximum deviation of %d',maxdeviation);
%         end 
        bond = bond.*[len len len];
    else            
        % [ bond ] = geometry2bond( cgrad, tgrad, tpgrad, len )        
        naa = length(c)+2;
        nfrags = naa-3;
        nx = size(c.dx,2);
        onx = ones(nx,1);        
        i123 = [1;2;3]; i132 = [1;3;2]; i231 = [2;3;1]; i321 = [3;2;1]; i312 = [3;1;2]; i213 = [2;1;3];
        
        if (verifyhandAD || useintlabAD)
            x = gradient([c;t;tp]);
            
            s = sqrt(1-x(1:naa-2).^2);
            
            bond = gradient(zeros(naa-1,3));
            bond(1,1) = 1;
            bond(2,1) = x(1);
            bond(2,2) = s(1); 
            ss = s(1:end-1).*s(2:end);     
        end
        % c(i) == x(i) (i \in 1..nfrags+1)
        % t(i) == x(nfrags+1+i) (i \in 1..nfrags)
        % tp(i) == x(2*nfrags+2+i) (i \in 1..nfrags)                                                
             
        % initialize x, get out of gradient object
        Ax = [c.x;t.x;tp.x];
        Dx = [c.dx;t.dx;tp.dx];
        
        % As = sqrt(1-c.^2)
        As1 = Ax(1:naa-2).^2; % c.^2
        Ds1 = 2*Ax(1:naa-2,onx).*Dx(1:naa-2,:);
        As2 = 1 - As1; % 1-c.^2
        Ds2 = -Ds1;
        As = sqrt(As2); % sqrt(1-c.^2)
        Ds = Ds2./(2*As(:,onx));
        
        % Ass = s_i * s_i-1
        Ass = As(1:end-1) .* As(2:end);
        Dss = Ds(1:end-1,:) .* As(2:end,onx) + As(1:end-1,onx) .* Ds(2:end,:);
        
        Abond = zeros(naa-1,3);
        Abond(1,1) = 1;
        Abond(2,1) = Ax(1); % c_1
        Abond(2,2) = As(1); % s_1 % TODO: we lack *sign(c(1))
        
        Dbond = zeros(naa-1,3,nx);
        Dbond(2,1,:) = Dx(1,:);
        Dbond(2,2,:) = Ds(1,:); % TODO: we lack *sign(c(1))

        for i=1:nfrags,            
            xci = i; % index of c_i in x
            xti = nfrags+1+i; % index of t_i in x
            xtpi = 2*nfrags+1+i; % index of tp_i in x
            
%             if (nargin == 5 && doAD && (x(xti).x ~= t(i) || x(xtpi).x ~= tp(i) || x(xci).x ~= c(i)))
%                 error('error');
%             end
           
            % Assemble our matrix A
            AA = [Abond(i:i+1,:); Abond(i,i231).*Abond(i+1,i312)-Abond(i,i312).*Abond(i+1,i231)];            
            DA = [Dbond(i:i+1,:,:); Dbond(i,i231,:).*Abond(i+1,i312,onx) + Abond(i,i231,onx).*Dbond(i+1,i312,:) - Abond(i,i312,onx).*Dbond(i+1,i231,:) - Dbond(i,i312,:).*Abond(i+1,i231,onx)];
                        
            % invert A: calculate determinant
            AdetA1 = AA([i123;i132],1).*AA([i231;i321],2);
            DdetA1 = AA([i123;i132],1,onx).*DA([i231;i321],2,:) + DA([i123;i132],1,:).*AA([i231;i321],2,onx);
            
            AdetA2 = [AA(i312,3);-AA(i213,3)];
            DdetA2 = [DA(i312,3,:);-DA(i213,3,:)];
            
            AdetA3 = AdetA1.*AdetA2;
            DdetA3 = DdetA1.*AdetA2(:,1,onx) + AdetA1(:,1,onx).*DdetA2;
            
            AdetA = sum(AdetA3);
            DdetA = sum(DdetA3,1);
            
            % invert A: calculate matrix for inversion. det(A)*Ainv == A
            AAinv = AA([2 3 3 1 1 2],i231).*AA([3 2 1 3 2 1],i312);
            DAinv = AA([2 3 3 1 1 2],i231,onx).*DA([3 2 1 3 2 1],i312,:) + DA([2 3 3 1 1 2],i231,:).*AA([3 2 1 3 2 1],i312,onx);
            
            AAinv = AAinv([1 3 5],:) - AAinv([2 4 6],:);
            DAinv = DAinv([1 3 5],:,:) - DAinv([2 4 6],:,:);
            
            % invert A: calculate Ainv = A/det(A)
            AAinv2 = AAinv/AdetA;
            DAinv2 = (DAinv*AdetA - AAinv(:,:,onx).*DdetA([1 1 1],[1 1 1],:))/AdetA^2;
            
            % form our vector b
            Ab = [Ax(xci)*Ax(xci+1) - Ax(xti)*Ass(i)
                  Ax(xci+1)
                  Ax(xtpi)*Ass(i)]';
            Db = zeros(nx,3);
            Db(:,1) = (Dx(xci,:).*Ax(xci+1,onx) + Ax(xci,onx).*Dx(xci+1,:) - Dx(xti,:).*Ass(i,onx) - Ax(xti,onx).*Dss(i,:))';
            Db(:,2) = Dx(xci+1,:)';
            Db(:,3) = (Dx(xtpi,:).*Ass(i,onx) + Ax(xtpi,onx).*Dss(i,:))';
            
            % multiply Ainv*b = bi == bond(i+2)
            Abi = zeros(1,3);
            Abi(1,1) = Ab*AAinv2(:,1);
            Abi(1,2) = Ab*AAinv2(:,2);
            Abi(1,3) = Ab*AAinv2(:,3);
            
            Dbi = zeros(nx,3);
            Dbi(:,1) = Db*AAinv2(:,1) + (Ab*squeeze(DAinv2(:,1,:)))';
            Dbi(:,2) = Db*AAinv2(:,2) + (Ab*squeeze(DAinv2(:,2,:)))';
            Dbi(:,3) = Db*AAinv2(:,3) + (Ab*squeeze(DAinv2(:,3,:)))';
                                    
            if verifyhandAD || useintlabAD,                 
                A = [bond(i:i+1,:); bond(i,i231).*bond(i+1,i312)-bond(i,i312).*bond(i+1,i231)];    
                detA = sum(A([i123;i132],1).*A([i231;i321],2).*[A(i312,3);-A(i213,3)]);

                Ainv = A([2 3 3 1 1 2],i231).*A([3 2 1 3 2 1],i312);
                Ainv = Ainv';
                Ainv = Ainv(:,[1 3 5]) - Ainv(:,[2 4 6]);
                Ainv = Ainv/detA;
                
                ssi = ss(i);
                xxcip = x(xci+1);
                b = [x(xci)*xxcip - x(xti)*ssi;xxcip;x(xtpi)*ssi];
                
                bi = (Ainv*b)';               
                
                graddiff = squeeze(bi.dx) - Dbi';
                if (norm(Abi-bi.x) > 1e-10) || any(graddiff(:) > 1e-10)
                    error('something went wrong');
                end
                
                nbi1 = bi.^2;
                nbi2 = sum(nbi1);                
                nbi = sqrt(nbi2);
                bond(i+2,:) = bi/nbi;                
            end
                                                
            % stabilize, divide bi by its norm            
            Anbi1 = Abi.*Abi;
            Dnbi1 = 2*Abi(onx,:).*Dbi;
            
            Anbi2 = sum(Anbi1);
            Dnbi2 = sum(Dnbi1,2);
            
            Anbi = sqrt(Anbi2);
            Dnbi = Dnbi2./(2*Anbi2);
                        
            Abond(i+2,:) = Abi'/Anbi;
            Dbond(i+2,:,:) = (Dbi*Anbi - Abi(onx,:).*Dnbi(:,[1 1 1]))'/Anbi^2;
                       
        end % for
        
        % multiply bonds by length
        len3 = [len len len];
        Abond = Abond.*len3;
        Dbond = Dbond.*len3(:,:,onx);
        
        % finally, repackage into gradient object
        r.x = Abond;
        r.dx = reshape(Dbond,3*(naa-1),nx);
        bond2 = gradient(r);
        
        if (useintlabAD)
            bond = bond.*len3;
        end
        if (verifyhandAD)            
            xdiff = abs(bond.x - bond2.x);
            dxdiff =  abs(bond.dx-bond2.dx);
            if (any(xdiff(:) > 1e-10) || any(dxdiff(:) > 1e-10))
                    error('something went wrong');
            end
        end
        
        if (~useintlabAD)
            bond = bond2;
        end
    end
    	    
    
end

