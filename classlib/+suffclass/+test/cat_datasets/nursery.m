% NURSERY Parser for the nursery data set.
function [ data ] = nursery()
nattrs = 9;

data = zeros(1,nattrs);

fid=fopen('nursery.data');
l = 1;
while 1
    tline = fgetl(fid);
    rline = tline;
    if ~ischar(tline) || strcmp(strtok(tline,','),'') ,   break,   end
    
    for i = 1:nattrs,
        [token, rline] = strtok(rline, ',');
        data(l,i) = catval(i, token);
        if (data(l,i) == 0)
            dprintf('problem parsing line #%i:\n"[%s]"\ntoken %i=%s',l,tline,i,token);
        end
    end
    %disp(tline)
    l = l+1;
end
fclose(fid);
 
end

function cat = catval(index, token)
%    parents        usual, pretentious, great_pret
%    has_nurs       proper, less_proper, improper, critical, very_crit
%    form           complete, completed, incomplete, foster
%    children       1, 2, 3, more
%    housing        convenient, less_conv, critical
%    finance        convenient, inconv
%    social         non-prob, slightly_prob, problematic
%    health         recommended, priority, not_recom
    cat = 0;
    switch(index) 
        case 1 
            switch(token)
                case 'usual' 
                    cat = 1;
                case 'pretentious'
                    cat = 2;
                case 'great_pret'
                    cat = 3;
            end
        case 2
            switch(token)
                case 'proper' 
                    cat = 1;
                case 'less_proper'
                    cat = 2;
                case 'improper'
                    cat = 3;
                case 'critical'
                    cat = 4;
                case 'very_crit'
                    cat = 5;
            end   
        case 3
            switch(token)
                case 'complete' 
                    cat = 1;
                case 'completed'
                    cat = 2;
                case 'incomplete'
                    cat = 3;
                case 'foster'
                    cat = 4;
            end 
        case 4
            switch(token)
                case '1' 
                    cat = 1;
                case '2'
                    cat = 2;
                case '3'
                    cat = 3;
                case 'more'
                    cat = 4;
            end 
        case 5
            switch(token)
                case 'convenient' 
                    cat = 1;
                case 'less_conv'
                    cat = 2;
                case 'critical'
                    cat = 3;
            end 
        case 6
            switch(token)
                case 'convenient' 
                    cat = 1;
                case 'inconv'
                    cat = 2;
            end   
        case 7
            switch(token)
                case 'nonprob' 
                    cat = 1;
                case 'slightly_prob'
                    cat = 2;
                case 'problematic'
                    cat = 3;
            end      
        case 8
            switch(token)
                case 'recommended' 
                    cat = 1;
                case 'priority'
                    cat = 2;
                case 'not_recom'
                    cat = 3;
            end 
        case 9
            switch(token)
                case 'not_recom' 
                    cat = 1;
                case 'recommend'
                    cat = 2;
                case 'very_recom'
                    cat = 3;
                case 'priority'
                    cat = 4;
                case 'spec_prior'
                    cat = 5;                    
            end 
    end
end