classdef Drawable < handle
    %DRAWABLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function draw(obj,h)            
            vis.ribbon3d(obj.getcoords.coords,10);
        end
        function drawnice(obj,h)                            
            vis.ribbon3dnice(obj.getcoords.coords);
        end
    end
    
    methods (Abstract)
        coords = getcoords(obj);
    end
    
end

