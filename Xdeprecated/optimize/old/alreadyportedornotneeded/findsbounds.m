function [clower,cupper,cplower,cpupper,zboth,zlower,zupper] = findsbounds(data,newclasses)
% I think these bounds should be predetermined, I'm using the max and min
% of ci and ti for aa's of type si found in the data set
% this program does not allow +- infinity as bounds, although such numbers
% occur. if bound actually is +-inf, sets to 0

if ~newclasses
    nums = 27;
else
    nums = 4; %HEFG classes
end

clower = repmat(inf,nums,1);
cupper = repmat(-inf,nums,1);
cplower = repmat(inf,nums,1);
cpupper = repmat(-inf,nums,1);
zlower = repmat(inf,nums,1);
zupper = repmat(-inf,nums,1);

global alphas;

for p = 1:size(data,2)
    if mod(p,100 )== 1
        p
    end
    bond = double(data{p}.bond);
    %normalize bond
    %bond = bond./repmat(sqrt(sum(bond.^2,2)),1,3);
    bond = rotatematrix(bond);
    seq = data{p}.seq;
    
    %class sequences
    if newclasses == 0 
        s = seq2sclasses(bond);
    else
        %skip this, s is in bond2geo, anyway
    end
    
    %finds geometry of normalized bond (geo.z, geo.c, geo.cp)
    if newclasses == 0
        geom = bond2geometry(bond); 
        geo = geometry2geometryz(geom,alphas);
    else
        geo = bond2geo(bond);
        s = geo.tcl;
    end
    
    sizes = size(seq,1)-3;
    for a = 1:sizes
        %gi = gamma(a);
        si = s(a);
        if newclasses == 0
           c = geom.c(a); %these are all *100
           cp = geom.cp(a); %these are all *100
           z = geo.zo(a); %this is not times 100
        else
           c = geo.co(a);
           cp = geo.cpo(a); 
           z = geo.tqu(a);
        end
        if ((c < clower(si)) && (c ~= -inf))
            clower(si) = c;
        end
        if ((c > cupper(si)) && ( c~= inf))
            cupper(si) = c;
        end
        if ((cp < cplower(si)) && (cp ~= -inf))
            cplower(si) = cp;
        end
        if ((cp > cpupper(si)) && (cp ~= inf))
            cpupper(si) = cp;
        end
        if ((z < zlower(si))) %&& (z ~= -inf))
            zlower(si) = z;
        end
        if ((z > zupper(si)))% && (z ~= inf))
            zupper(si) = z;
        end
   
    end % end for a
end % end for p

%get rid of infinities
for i = 1:nums
    if abs(cupper(i)) == inf
        cupper(i) = 1;
    end
    if abs(clower(i)) == inf
        clower(i) = -1;
    end
    if abs(cpupper(i)) == inf
        cpupper(i) = 1;
    end
    if abs(cplower(i)) == inf
        cplower(i) = -1;
    end
    if zupper(i) == -inf %occurs if the class s doesn't exist (long 0 columns in confus)
        zupper(i) = 100;
    end
    if zlower(i) == inf  %occurs if the class s doesn't exist (long 0 columns in confus)
        zlower(i) = -100;
    end
end
    
zboth = zeros(nums,1);
%take absolute max as zboth
for i = 1:nums
    if abs(zlower(i)) >= abs(zupper(i))
        zboth(i) = abs(zlower(i));
    else
        zboth(i) = abs(zupper(i));
    end
end

%adjust c, cp bounds so that they are compatible
% although I really don't get why they're not compatible

for i = 2:nums
    if clower(i) >= cpupper(i-1)
        %now should I raise cpupper or lower clower?
        clower(i) = cplower(i-1);
    end
    if cupper(i) <= cplower(i-1)
        cupper(i) = cpupper(i-1);
    end
end

        