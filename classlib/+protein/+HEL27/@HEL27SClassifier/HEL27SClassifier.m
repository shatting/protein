classdef HEL27SClassifier < protein.SClassifier
    %HEL27CLASSIFIER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods   
        
    end
    
    methods
        
        function obj = HEL27SClassifier()
           [bounds, alphas] = class.HEL27SClassifier.bounds;
           obj = obj@class.SClassifier(bounds,alphas); 
        end
        
        cl = classify(obj,dataset)           
        
        function featurenames = requiredfeaturenames(obj)
            featurenames = {{'c','cp','beta'}};
        end  
        
        function ncl = getncl(obj)
            ncl = 27;
        end
    end
    methods (Static)
                
        function [bounds,alphas] = bounds
           
            % TODO: it might be that the c, cp classes given in alphas are not in
            % correspondence to the ordering in bounds
            % TODO: make, cbounds, cpbounds, betabounds

            % these are c and cp dependent bounds in the torsion angles
            bounds = [-40  60  inf
                      -70 -14  50
                      -88 -18  40
                      -80  -5  30
                      -70  30 inf
                      -80  20 inf
                      -70  30 inf
                      -50  40  75         
                      -60  20  78];

            % these are the c and cp dependent torsion class centers
            % TODO: check this. ERROR: these were meant to be the ccp dependent torsion class BOUNDARIES
            alphas = [[1 1 10 -10 inf];...
                      [1 2 -42 18 90];...
                      [1 3 -53 11 76];...
                      [2 1 -42.5 12.5 75];...
                      [2 2 -20 80 inf];...
                      [2 3 -30 70 inf];...
                      [3 1 -20 80 inf];...
                      [3 2 -5 57.5 -85.7];...
                      [3 3 -20 42 -91]]; %not sure about the 42, should make a function to calculate these, since they'll probably change            
 
        end
    end
    
end

