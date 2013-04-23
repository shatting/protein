classdef Classification < data.Dataset
    %CLASSIFICATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        classifier
        cl
        covmodel
    end
    
    methods
        
        % copies the required features
        function obj = Classification(dataset,inclassfier)            
            obj.classifier = inclassfier;
        
            % decide which feature set to use
            fnames = inclassfier.requiredfeaturenames;
            for i=1:length(fnames)
                if (dataset.haveallfeatures(fnames{i}))
                    fnames = fnames{i};
                    break;
                end
            end
            % TODO: ? only give the required featureset to classify routine
            obj.cl = inclassfier.get_classes(dataset);
            
            obj.initialized(dataset.getdata(fnames),fnames);        
        end
       
        
    end
    
    methods (Access=protected)
                
        % method stub so we can instantiate. initialization is done in
        % constructor
        function initialize(dataset)
        end
    end
    
end

