function [ t, tp ] = ttptransform_rational( z, sseq, cosbetabar, sinbetabar )
%TTPTRANSFORM_RATIONAL Summary of this function goes here
%   Detailed explanation goes here

sinbetabari = sinbetabar(sseq);
cosbetabari = cosbetabar(sseq);
                 
zsquared = z.^2;            
sinbetaminusbetabar = 2*z./(1+zsquared);
cosbetaminusbetabar = (1-zsquared)./(1+zsquared);

t = -sinbetabari.*sinbetaminusbetabar + cosbetabari.*cosbetaminusbetabar;
tp = cosbetabari.*sinbetaminusbetabar + sinbetabari.*cosbetaminusbetabar;

end

