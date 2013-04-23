function [ sc ] = stringcellremove( sc, idx )
%STRINGCELLREMOVE 

if (idx < 1 || idx > length(sc))
    stringcell = sc
    error('string cell index %i out of bounds.',idx);
end

if idx==1,
    sc = sc(2:end);
elseif idx == length(sc),
    sc = sc(1:end-1);
else
    sc = sc([1:idx-1,idx+1:end]);
end

end

