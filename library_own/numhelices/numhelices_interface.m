function [ numhelices ] = numhelices_interface( data, new_deal_data, sclasses )
%NUMHELICES_INTERFACE the interface of this package

[ numscln, midx ] = numsclasseslengthn( data, new_deal_data, sclasses );

numhelices = numhelices_filter(numscln);

end
