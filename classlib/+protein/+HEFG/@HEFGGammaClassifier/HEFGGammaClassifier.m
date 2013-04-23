classdef HEFGGammaClassifier < class.PotentialClassifier & protein.GammaClassifier
    %HEFGGAMMACLASSIFIER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hefg_trainingrun
    end
    
    methods 
        function  this = HEFGGammaClassifier(hefg_trainingrun)            
            this = this@class.PotentialClassifier(hefg_trainingrun.options.featfun, hefg_trainingrun.class_result.potential_values);
            this.feature_fn = hefg_trainingrun.options.featfun;
            this.hefg_trainingrun = hefg_trainingrun;
        end    

       function ncl = getncl(this)
          ncl = this.potential_values.potential.suff.ncl;          
       end

       function names = requiredfeaturenames(this)
            names = {{'aa1','aa2','aa3','aa4'}}; 
       end    
             
    end
    
    methods (Access=protected)

    end

end

