function [ threefrags ] = threefrags_interface( doclassify)%,rawdata )
%[ threefrags ] = threefrags_interface( doclassify )
% get desired potential for (aa1,aa2) by potc(pair2cl(aa1,aa2))    

% compare with stephan's paper!!! 11 of the entries in the matrix are used
% as features - make sure to get them in the correct order!! the potential
% is in potc - L is not R^-T, LL^T = inv(cov) or something
% everything has to be amplized 

%getrawdata

% we look at pair distances 2:0.5:7
load rawdata_r_2_to_7;

% get suffdata struct with aa pair suffs, weighted by distance function
raw2suff_w;

if doclassify,
    classify_w
    cl = cl(:,end);
else    
    cl=1:400;
end

ij = sys10toX(1:400,[20,20]);
threefrags.pair2cl = full(sparse(ij(:,1),ij(:,2),cl,20,20));
threefrags.potc = suffs2pots(suffdata,cl);


% these are used to calculate the weights
wtfn_norm
wtfn
maxwt
disp('rawdata.rbins = ')
disp(rawdata.rbins);

disp('th.potc(1).cov*th.potc(1).L''*th.potc(1).L = I');

end
