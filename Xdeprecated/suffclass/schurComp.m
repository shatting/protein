function [ CS ] = schurComp( C, idx )
%SCHURCOMP returns the schur complement of C, denoted by C_[idx]

idx2 = setdiff(1:size(C,1),idx);

CS = C(idx,idx) - C(idx,idx2)*inv(C(idx2,idx2))*C(idx2,idx);

end
