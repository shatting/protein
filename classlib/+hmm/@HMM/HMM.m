classdef HMM
    %HMM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sclassifier
        gclassifier        
        sgconfus
    end
    
    methods
        function obj = HMM(sclassifier, gclassifier, dataset)
            obj.sgconfus = hmm.getsgconfus(dataset,sclassifier,gclassifier);
            obj.sclassifier = sclassifier;
            obj.gclassifier = gclassifier;
        end
        
        function hmmseqs = getnbest(obj,chain,nfound)
            gammaseq = obj.gclassifier.classify(chain);
            hmmseqs = hmm.getnbest(gammaseq,obj.sgconfus,nfound);
        end
        
    end
    
end

