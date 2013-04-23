function [ wmeankld, meankld ] = kldivqualitymeasure( target, cl, verbose )
%KLDIVQUALITYMEASURE Summary of this function goes here
%   Detailed explanation goes here

conf = suffclass.utils.getfreq([target cl]);

    if (verbose)
        pairkld = protein.kldiv.kldivpairs(conf);

        svdpair = svd(pairkld)'
        svratio = max(svdpair)/min(svdpair(svdpair>0))
        svdpairsym = svd(.5*pairkld*pairkld')'
        svsymratio = max(svdpairsym)/min(svdpairsym(svdpairsym>0))

        pr = sum(conf)/sum(sum(conf));
        kld = protein.kldiv.kldiv(conf,sum(conf,2))
        meankld = mean(kld)
        wmeankld = sum(kld.*pr)
    else
        pr = sum(conf)/sum(sum(conf));
        kld = protein.kldiv.kldiv(conf,sum(conf,2));
        meankld = mean(kld);
        wmeankld = sum(kld.*pr);
    end

end
