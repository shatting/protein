classdef Geometry < data.Dataset & vis.Drawable
    %GEOMETRY the primary representation for a chain geometry. everything
    %here is rotational invariant.
    %
    % needed: 
    %   - convert from/to:  * original bond data
    %                       * optimization results (z)
    %                       * normalized bonds (we need to know that those
    %                           are normalized)
    %                       * other geometry object
    %   - methods to:       * get s class sequence
    %                       * rotate
    %                       * shift
    %                       * normalize
    %                       * center
    %                       * register
    
    
    properties (Constant)
        version = 1;
    end
    
    properties (SetAccess = protected, GetAccess=protected)        
    end
    
    % these are fields initialized 
    properties (SetAccess=protected,Transient)        
        c;
        s;
        t;
        tp;     
        svdbb;
        naa;
        nfrags;     
        beta;    
        len;
    end
        
    methods (Access = private)
        % initializer method from bond
        [ obj ] = fromdatabond(obj, bond)
    end
      
    methods
      % constructor
      function obj = Geometry(bond)                 
        obj = fromdatabond(obj,bond);
      end      
      
      % Drawable implementation
      function [coords, Acond] = getcoords(obj)
            if (nargout > 1)
                [bond, Acond] = obj.getbond;
            else
                bond = obj.getbond;
            end
            coords = geom.Coords(bond);
      end           
      
    end
    methods (Access = protected)
        % Dataset implementations
        % initialize(obj)
        % extract fields from geomstruct and returns a matrix containing (all/idx) fragments' 
        % features 
        % INPUT:    geomstruct   
        %           fields       cell array of fieldnames, all fields if omitted or
        %                        empty. special treatment of cp/sp: if c is in
        %                        fieldnames, cp will be added although it is no
        %                        fieldname. same applies to s. this is due to the
        %                        fact that c and s dont directly correspond to
        %                        fragments (cp and sp need to be shifted by 1)
        %           idx          fragment indices
        % OUTPUT:
        %           fielddata    (nfrags or |idx| x nfields) matrix
        %           fields       fields actually used. this is useful if no fields
        %                        were given
        % uses: (none)
        function initialize(obj)            
                  
            % get field names
%             mc = ?geom.Geometry;
%             names = cell(0);
%             for i=1:length(mc.Properties)
%                 if (~mc.Properties{i}.Constant)
%                     names = [names,mc.Properties{i}.Name];
%                 end
%             end
%             % remove those names not suitable for feature creation
%             names = setdiff(names,{'bondlen','res','r','bond','svdbb','naa','nfrags','data','featurenames','len','ndata','nfeats'});
%             names = [names,'cp','sp'];
            %safer: 
            names = {'c','cp','s','sp','t','tp','beta'};
            names = sort(names);
            
            nfeats = length(names);
            idx = 1:length(obj.beta);
            features = zeros(length(idx),nfeats);

            for i=1:nfeats,
                if strcmp(names{i},'cp'),
                    temp = obj.c(2:end);
                elseif strcmp(names{i},'sp'),
                    temp = obj.s(2:end);
                else
                    temp = obj.(names{i});
                end
                features(:,i) = temp(idx);
            end

            obj.initialized(features,names);
        end
        
        function [ bond, Acond ] = getbond( obj )
            if (nargout > 1)
                [bond, Acond] = geom.geometry2bond(obj.c,obj.t,obj.tp,obj.len);
            else
                bond = geom.geometry2bond(obj.c,obj.t,obj.tp,obj.len);
            end
        end
      end
      
end    