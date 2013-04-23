classdef TrainingRun < handle
    % suffclass.TrainingRun
    %   abstract base class for all training runs
    %
    % PROPERTIES
    %   options
    % 
    
    properties
        options    
    end
    
    methods (Abstract)
        run(this)
    end
    
    methods
                
    end
    
end

