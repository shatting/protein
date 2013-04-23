classdef HEL27GammaClassifier < protein.GammaClassifier & class.PotentialClassifier
    %HEL27GAMMACLASSIFIER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
       function  obj = HEL27GammaClassifier(potential)
            obj = obj@class.PotentialClassifier(potential);
       end
                
       function ncl = getncl(obj)
          ncl = max(obj.potential.final.cl(:,end));
       end
       
       function names = requiredfeaturenames(obj)
            names = {{'aa1','aa2','aa3','aa4'}}; 
       end
    end
            
    methods (Access=protected)
        cl = getclasses(obj,dataset)
    end
        
end

