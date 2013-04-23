classdef crossval
    %CROSSVAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nfold
        perm
        foldidxstart
        n
    end
    
    methods
        function obj = crossval(n,nfold,lastsmaller)
           obj.nfold = nfold; 
           obj.perm = randperm(n);                      
           testsize = round(n/nfold);           
           obj.foldidxstart = 1:testsize:n;
           if (length(obj.foldidxstart) ~= nfold + 1)
               obj.foldidxstart(end+1) = n + 1;
           end
           obj.foldidxstart(end) = n+1;
           obj.n = n;
        end
        
        function [testidxset, trainidxset] = foldidxset(obj, foldnum)
            if foldnum>obj.nfold, error('CV: Fold number is greater than number of Folds'); end
            testpermidx = obj.foldidxstart(foldnum):obj.foldidxstart(foldnum+1)-1;
            trainpermidx = setdiff(1:obj.n,testpermidx);
            testidxset = obj.perm(testpermidx);            
            trainidxset = obj.perm(trainpermidx);            
        end
        
        function display(obj)
            dprintf('%i-fold cross validation, test set sizes:',obj.nfold);
            suffclass.utils.tightmat(diff(obj.foldidxstart));
        end
    end
    
end

