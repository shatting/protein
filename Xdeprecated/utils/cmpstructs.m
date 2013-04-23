function [ err ] = cmpstructs( str1, str2 )
%CMPSTRUCTS Summary of this function goes here
%   Detailed explanation goes here

err = 0;
f = fields(str1);
for i=1:length(f),
    field1 = eval(['str1.',f{i}]);
    field2 = eval(['str2.',f{i}]);
    
    dfield = field1-field2;
    errf = sum(abs(dfield(:)));
    
    dprintf('error in field %s: %d', f{i}, errf);
    
    err = err + errf;
end

end
