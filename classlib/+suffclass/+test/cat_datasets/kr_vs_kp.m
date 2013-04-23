%KR_VS_KP Parser for kr_vs_kp problem set.

function [ data ] = kr_vs_kp()

nattrs = 37;

data = zeros(1,nattrs);

fid=fopen('kr-vs-kp.data');
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

for i=1:36,
   u = unique(data(:,i))';
   m = max(u);
   u = [u setdiff(1:m,u)];
   [s, p] = sort(u);
   data(:,i) = p(data(:,i));
   dprintf('i: ');disp(unique(data(:,i)));
end
 
    function cat = catval(index, token)
        cat = 0;
        if (index>=1 && index <=36) 
                switch(token)
                    case 'f' 
                        cat = 1;
                    case 't'
                        cat = 2;
                    case 'l'
                        cat = 3;
                    case 'n'
                        cat = 4;
                    case 'w'
                        cat = 5;
                    case 'b'
                        cat = 6;
                    case 'g'
                        cat = 7;
                    otherwise
                        dprintf('unknown token:%s',token);
                end
        end
        if index == 37,
                switch(token)
                    case 'won' 
                        cat = 1;
                    case 'nowin'
                        cat = 2;
                end 
        end
    end

end