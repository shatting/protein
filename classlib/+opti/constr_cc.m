function [ cc, distances, cci ] = constr_cc( coords )
%CONSTR_CC Summary of this function goes here
%   Detailed explanation goes here
    % we want to calculate in angstöm
    coords = coords/10;
    n = size(coords,1);
    testhandgradient = 0;
    dograd = 0;
    
    if 0,
        dists = zeros(n,n);
        mask = triu(dists+1,1);
        for i=1:n,
            for k=i+1:n,
                dvec = coords(i,:) - coords(k,:);
                dists(i,k) = sum(dvec.^2);
            end
        end
        dsquared = dists(logical(mask(:)));   
        cc = sum(1./(1e-4 + max(dsquared,1e-4)));
    else
        if (isa(coords,'gradient'))
            if (testhandgradient)
                dvecs = gradient(zeros(n*(n-1)/2,3));   % distance vectors                  
            end
            nxgrad = size(coords.dx,3);
            onxgrad = ones(1,nxgrad);
            Advecs = zeros(n*(n-1)/2,3);
            Ddvecs = zeros(n*(n-1)/2,3,nxgrad);
            Acoords = coords.x;
            Dcoords = coords.dx;
            dograd = 1;
        else
            dvecs = zeros(n*(n-1)/2,3);   % distance vectors                    
        end
        if (isa(coords,'opti.Coords'))
            coords=coords.coords;
        end
        
        idxs = 1;
        for i=1:n-1,            
            idxe = idxs + n-i -1;
            on = i*ones(n-i,1);
            idxminus = i+1:n;
            if (testhandgradient || ~dograd)
                dvecs(idxs:idxe,:) = coords(on,:) - coords(idxminus,:);
            end
            if (dograd)
                Advecs(idxs:idxe,:) = Acoords(on,:) - Acoords(idxminus,:);
                Ddvecs(idxs:idxe,:,:) = Dcoords(on,:,:) - Dcoords(idxminus,:,:);
            end
            idxs = idxe +1;
        end        
        
        if (testhandgradient || ~dograd)
            dsquared = sum(dvecs.^2,2);  
            cci = 1./(1e-4 + dsquared);
            cc = sum(cci);
        end
        
        if (dograd)
            Adsquared1 = Advecs.^2;
            Ddsquared1 = 2*Ddvecs.*Advecs(:,:,onxgrad);
            
            Adsquared2 = sum(Adsquared1,2) + 1e-4;
            Ddsquared2 = permute(sum(Ddsquared1,2),[1 3 2]);                        
            
            Acci = 1./Adsquared2;
            Acci2 = Acci.^2;
            Dcci = -Ddsquared2.*Acci2(:,onxgrad);
            cc2.x = sum(Acci);
            cc2.dx = sum(Dcci);
            cc3 = gradient(cc2);
        end
                
        if (testhandgradient && dograd)
            dd = cc - cc3;
            if max(abs(dd.x)) > 1e-8 || max(abs(dd.dx)) > 1e-10,
                error('');
            end            
        end
        
        if (dograd)
            cc = cc3;
        end
        
    end
    if (nargout > 1)        
        distances = sqrt(dsquared);
    end
end

