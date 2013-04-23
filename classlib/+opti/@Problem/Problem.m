classdef Problem
    %PROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        priseq
        sclassifier
        gclassifier
    end
    
    methods
        function obj = problem(priseq, sclassifier, gclassifier)
            obj.priseq = priseq;
            obj.sclassifier = sclassifier;
            obj.gclassifier = gclassifier;
        end
        
        
        
    end
    
end

