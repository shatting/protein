function [ features ] = fragdata_getfeaturesbyname( fragdata, featurenames )
%[ features ] = featuresfromfragdata( fragdata, featurenames )
% extract features with names featurenames from fragdata

idxf = zeros(1,length(featurenames));

for i=1:length(featurenames),
    idxf(i) = stringcellindexof(fragdata.featurenames,featurenames(i));
end

features = fragdata.data(:,idxf);

end

