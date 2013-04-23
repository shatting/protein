classdef SimpleDataset < data.Dataset
    % data.SimpleDataset < data.Dataset
    %   Most basic dataset.
    
    properties
    end
    
    methods
        function obj = SimpleDataset(data,featurenames)
            obj.initialized(data,featurenames);
        end

    end
    
    methods (Access=protected)        
        function initialize(dataset)        
            
        end
    end
    
end

