function [ names ] = prednames( conf )
%NEW_DEAL_HEFG_TAGETNAMES Summary of this function goes here
%   Detailed explanation goes here

full = 0.9;
prime = 0.5;

[nt np] = size(conf);

conf_max = conf./repmat(max(conf),nt,1);

names = cell(np,1);
maxname = 0;

for i=1:np,
    siz = conf_max(:,i);
    rsiz = round(siz*100);  
    names{i} = [compname('H',siz(1)),compname('E',siz(2)),compname('F',siz(3)),compname('G',siz(4))];
    maxname = max(maxname,length(names{i}));
end

for i=1:np,   
    siz = conf_max(:,i);
    rsiz = round(siz*100); 
    pad = ' ';
    pad = repmat(pad,1,maxname-length(names{i})+1);
    names{i} = [names{i} ,sprintf('%s[%i,%i,%i,%i]',pad,rsiz(1),rsiz(2),rsiz(3),rsiz(4))]; 
end

    function n = compname(name,val)
        if (val>=full)
            n = name;
        elseif val>=prime
            n = [name,''''];
        else
            n='';
        end
    end

end
