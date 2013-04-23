classdef Chain < geom.Coords & data.Dataset
    % geom.Chain < geom.Coords & data.Dataset
    %   
    
    properties (SetAccess = private)
        name;
        seq;        
        res;
    end
    
    methods
        function this = Chain(chainstruct)           
            this = this@geom.Coords(chainstruct.bond);
            this.seq = chainstruct.seq;
            this.name = chainstruct.name;
            this.res = chainstruct.res;
        end
                
        function display(this)
            display@data.Dataset(this);
            if (length(this) == 1)
                dprintf('\nChain "%s", %i AA, res: %.2fAng', this.name, length(this.seq), this.res);
            end
        end
    end
    
    methods (Access = protected)
        
        % geom.Dataset implementations
        function initialize(dataset)            
            featurenames = [dataset.geometry.getfeaturenames(),'aa1','aa2','aa3','aa4'];
            
            % TODO: problem: data duplication
            odata = [dataset.geometry.getdata() seq2fragseq(dataset.seq)];        
        
            dataset.initialized(odata,featurenames);
        end         
    end
    
end

