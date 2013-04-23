function [ prot_s, prot_gamma, pot ] = ginny_interface( RTdata, geom, class_info )
%GINNY_INTERFACE 
% >>load fragdb
% >>load testclassify
% >>load RT071127
% >>geom = fragdb2linear(fragdb);
% >>[ s, gamma, pot ] = ginny_interface( data, geom, class_info )

prot_s = seq2sclasses(RTdata{1}.seq,RTdata{1}.bond); % protein s sequence
prot_gamma = seq2gammaclasses(RTdata{1}.seq,class_info); % protein gamma sequence

all_s = HEL27_anal08(geom);
all_gamma = class_info.cl(:,end);

pot = opticov( geom(:,1:3), all_s, all_gamma );

end
