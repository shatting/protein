classdef SClassification < class.Classification
    %SCLASSIFICATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function obj = SClassification(dataset,inclassfier)            
            obj = obj@class.Classification(dataset,inclassfier);                               
        end
        
    end
    
end

