classdef HEFGSClassifier < protein.SClassifier
    %HEFGSCLASSIFIER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        diagbds
    end
    
    methods
        function this = HEFGSClassifier(diagbds)            
            this = this@protein.SClassifier([]);                   
            this.diagbds = diagbds;
        end
                
        
        function names = requiredfeaturenames(this)
            names = {{'beta'},{'t','tp'}};
        end
        
        function ncl = getncl(this)
            ncl = 4;
        end
        
        function [sane, errors] = test(obj)
            dprintf('- HEFGClassifier(diag=%i).m',obj.diagbds);
            ds = data.FeatureDB.db.getdatasetclone();
            ds.removefeature('beta');
            clttp = class.Classification(ds,obj);
            ds = data.FeatureDB.db.getdatasetclone();
            ds.removefeature('t');
            ds.removefeature('tp');
            clbeta = class.Classification(ds,obj);
            errors = find(clbeta.cl~=clttp.cl);
            sane = sum(errors) == 0; 
            if sane,    
                dprintf('PASSED.');
            else
                dprintf('PROBLEMS FOUND!');
            end
        end
        
        % tc = gettorsioncenters(obj)
        % returns in tc(i) the center of class i beta angles for the z
        % transform
        function tc = gettorsioncenters(obj)
            if (obj.diagbds)
                tc = [pi/2,pi,-pi/2,0]; 
                % have to be careful transforming in E class
                % beta' = 2*pi+beta if beta < 0, then centering around pi makes sense
            else
                tc = [pi/4,-3*pi/4,3*pi/4,-pi/4];
            end
        end
        
        % H=1, E=2, F=3, G=4        
        function distm = getdistm(obj)
            if (obj.diagbds)
                distm = [0  1  2  1
                         1  0  1  2
                         2  1  0  1
                         1  2  1  0];
            else
                distm = [0  2  1  1
                         2  0  1  1
                         1  1  0  2
                         1  1  2  0];
            end
        end
        
        function showstuff(this)

            ttp = data.FeatureDB.db.getdata({'t','tp'});
            c = data.FeatureDB.db.getdata({'c'});
            cp = data.FeatureDB.db.getdata({'cp'});

            sclassifier = class.HEFGSClassifier(this.diagbds);
            hefg = sclassifier.classify(data.FeatureDB.db);

            radii = rand(size(ttp,1),1);
            radii = sqrt(repmat(radii,1,2));

            ttp_scaled = ttp.*radii;

            markersize = 0.01;
            colmap = colormap(class_colormap(4));
            subx = 2;
            suby = 2;

            figure;
            %subplot(suby,subx,[1 2 3 4]);
            plot(ttp_scaled(hefg==1,1),ttp_scaled(hefg==1,2),'.','MarkerSize',markersize,'Color',colmap(1,:))
            hold on;
            plot(ttp_scaled(hefg==2,1),ttp_scaled(hefg==2,2),'.','MarkerSize',markersize,'Color',colmap(2,:))
            plot(ttp_scaled(hefg==3,1),ttp_scaled(hefg==3,2),'.','MarkerSize',markersize,'Color',colmap(3,:))
            plot(ttp_scaled(hefg==4,1),ttp_scaled(hefg==4,2),'.','MarkerSize',markersize,'Color',colmap(4,:))
            %title('\beta, radius radomized');
            axis off;
            axis equal;
            hold off;

            figure
            subplot(suby,subx,[1]);
            vis.class_barhist( data.FeatureDB.db.getdata({'beta'}), hefg, 50, 0, 4, 0, 1, {'H','E','F','G'}, 0);
            title('Histogram of \beta');

            subplot(suby,subx,[2]);
            %class_barhist( hefg, hefg, 4, 0, 4, 0, 1, {'H','E','F','G'}, 0);
            vis.class_barhist(hefg,hefg,1:4,0,4,0,1,[], 0);
            title('s-class Frequencies');

            subplot(suby,subx,3);
            vis.class_barhist( c, hefg, 50, 0, 4, 0, 1, [], 0);
            title('Histogram of cos(\alpha)');

            subplot(suby,subx,4);
            vis.class_barhist( c, hefg, 50, 0, 4, 1, 1, [], 0);
            title('Scaled Histogram of cos(\alpha)');
            axis equal
        end
    end
    
    methods (Access=protected)
        cl = getclasses(obj,dataset)
    end
    
end

