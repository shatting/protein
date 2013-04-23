classdef Coords < handle & vis.Drawable
    % +geom@Coords
    %   Encapsulates one geometric configuration
    
    properties (SetAccess = protected, GetAccess = public)
        geometry
        bond   %
        x0     %
        coords % primary data store
        normalizerotm
    end
    
    properties (Transient, SetAccess = protected, GetAccess = protected)
        
    end   
    
    methods (Access = protected)
        function isotransformed(obj,bond,x0) % rotated, shifted
            obj.bond = bond;     
            obj.x0 = x0;
            obj.coords = geom.bond2coords(obj.bond,x0);                        
        end        
        
        function transformed(obj,bond,x0) % lengths may have changed too
            obj.isotransformed(bond,x0);            
            obj.geometry = geom.Geometry(bond);
        end        
    end
    
    methods
        function obj = Coords(bond,x0)
            bond = double(bond);
            if (nargin<2)
                x0 = [0 0 0];
            end
            obj.transformed(bond,x0);            
        end
               
        % Drawable implementations
        function coords = getcoords(obj)
            coords = obj;
        end
                
        function b = isunitlength(obj)
            b = abs(1 - max(sqrt(sum(obj.bond.^2,2)))) < 1e-10;
        end
        
        function b = isnormalizedrotation(obj)
            b = (norm(obj.bond(1,2:3)) < 1e-8) && obj.bond(2,3) < 1e-8;
        end
         
        function b = iscentered(obj)
            b = norm(mean(obj.coords,1)) < 1e-8;
        end
        
        function rotm = normalizerotation(obj)
            %[ rotm ] = bond_rotate( bond ) rotate bonds to canonical
            %configuration.
                        
            [bondr, rotm] = geom.bond_normalizerotation(obj.bond);                    
            obj.normalizerotm = rotm;
            obj.isotransformed(bondr,obj.x0);
        end
               
        function ormsd = rmsdto(obj,other)                                  
            ormsd = geom.rmsd(obj.coords,other.coords);
        end
        
        function [ormsd,Q,q] = minrmsdto(obj,other)                                  
            [Q,q,ormsd] = geom.coords_register(other.coords,obj.coords);
        end
        
        function rmsds = progressivermsdto(obj,other)
            rmsds = geom.rmsd_progressive(obj.coords,other.coords);
        end
        
        function registerreportedrmsd = registerto(obj,other)
            [Q,q,registerreportedrmsd] = geom.coords_register(other.coords,obj.coords);
            
            obj.rotateby(Q');
            obj.shiftby(q);
        end
                
        % iso transforms
        function shiftby(obj,vector)            
            obj.isotransformed(obj.bond,obj.x0 + vector);
        end
        
        function centermean(obj)
            obj.isotransformed(obj.bond,-mean(obj.coords,1));
        end
        
        function setx0to(obj,x0)
            obj.isotransformed(obj.bond,x0);
        end
        
        function rotateby(obj,Q)            
            obj.isotransformed(obj.bond*Q,obj.x0);
        end                             
        
    end
    
    methods (Static)
        
        %%% ----------------- TEST ROUTINES -------------------------------
        %%% ---------------------------------------------------------------
        
        function test_viewribbons(chnum)
            ch = data.GeomDB.ch{chnum};
            figure;
            %subplot(2,2,1);
            title('original');
            ch.drawnice;
            %subplot(2,2,2);            
            figure;
            title('normalized');
            ch.normalizerotation;
            ch.drawnice;
            %subplot(2,2,3);            
            figure;
            title('from Geometry');
            ch.geometry.drawnice;
        end
        
        function [problems,rmsds] = test_registerto(chnum,detailed)
            ch = data.GeomDB.db.chains{chnum};
            if (ch.isnormalizedrotation)
                data.GeomDB.reloadch(chnum);
                ch = data.GeomDB.db.chains{chnum};
            end
            newch = geom.Chain(ch);
            
            newch.normalizerotation;
            newch.registerto(ch);
            
            rmsds = ch.progressivermsdto(newch);
            
            problems = max(rmsds) > 1e-8;
            
            if ((nargin > 1 && detailed))                
                dprintf('- .registerto test. chnum = %i, name=%s, naa = %i, res = %.2f --',chnum,ch.name,length(ch.seq),ch.res);                                               
                dprintf('* rmsd of original to registerd = %g (should be < 1e-5)',rmsds(end));
                dprintf('* max progressive rmsd of original to registerd = %g (should be < 1e-5)',max(rmsds));
                
                plot(rmsds/max(rmsds),'+r');
                title(sprintf('.registerto test, chain %i',chnum));
                legend2 = sprintf('progressive rmsds (max=%g)',max(rmsds)); 
                legend({legend2});                
            end
        end
        
        function [problems,rmsds,registerreportedrmsd] = test_normalize(chnum,detailed)
            ch = data.GeomDB.db.chains{chnum};
            if (ch.isnormalizedrotation)
                data.GeomDB.reloadch(chnum);
                ch = data.GeomDB.db.chains{chnum};
            end
            
            coordsbefore = ch.getcoords;
            rotm = ch.normalizerotation;            
            coordsafter = ch.getcoords;
            
            % check if it is still the same chain
            [Q,q,registerreportedrmsd] = geom.coords_register(coordsbefore,coordsafter);
            
            registeredcoordsafter = geom.Coords(coordsafter.bond);
            registeredcoordsafter.rotateby(Q');
            registeredcoordsafter.shiftby(-q);            
            
            rmsds = geom.rmsd_progressive(coordsbefore,registeredcoordsafter);
            
            problems = max(rmsds)>1e-8 || ~ch.isnormalizedrotation;
            
            if ((nargin > 1 && detailed))                
                dprintf('- normalizerotation test. chnum = %i, name=%s, naa = %i, res = %.2f --',chnum,ch.name,length(ch.seq),ch.res);                               
                dprintf('* register: rotmcondition = %g',cond(Q));
                dprintf('* register: reportedrmsd = %g',registerreportedrmsd);                
                dprintf('* rmsd of original to registerd = %g (should be < 1e-5)',rmsds(end));
                dprintf('* max progressive rmsd of original to registerd = %g (should be < 1e-5)',max(rmsds));
                                    
                plot(rmsds/max(rmsds),'+r');
                title(sprintf('normalizerotation test, chain %i',chnum));
                legend2 = sprintf('progressive rmsds (max=%g)',max(rmsds)); 
                legend({legend2});                
            end
        end
        
        function [problems, rmsds, maxAcond] = test_geomreconstruction(chnum,detailed)            
            ch = data.GeomDB.db.chains{chnum};            
            if (ch.isnormalizedrotation)
                data.GeomDB.reloadch(chnum);
                ch = data.GeomDB.db.chains{chnum};
            end
            
            ch.normalizerotation;            
            ch.setx0to([0 0 0]);
            
            % reconstruction is implicit normalization
            % TODO: half of the chains get rotated configuration around 
            % x axis. reconstruction and normalization do not recognize
            % this, but registration does. maybe do normalization through registering to
            % [[1 0 0];[1 1 0];[0 0 0]] or something similar?
            [reccoords, Acond] = ch.geometry.getcoords;
            reccoords.registerto(ch);
            rmsds = ch.progressivermsdto(reccoords);
            
            problems =   rmsds(end) > 1e-8 || max(Acond) > 50;                        
                       % max(abs(lengthdeviations)) > 1e-10 ||
             
            maxAcond = max(Acond);
            
            if ((nargin > 1 && detailed))                
                dprintf('- normalizerotation+reconstruction test. chnum = %i, name=%s, naa = %i, res = %.2f --',chnum,ch.name,length(ch.seq),ch.res);                                
                %dprintf('* reconstruction: max length deviation = %g',max(abs(lengthdeviations)));
                dprintf('* reconstruction: max matrix condition = %g',max(Acond));                
                dprintf('* rmsd of rotated to reconstructed = %g (should be < 1e-8)',rmsds(end));                
                                
                plot(rmsds/max(rmsds),'+r');
                hold on;
                %plot(lengthdeviations/max(lengthdeviations));
                plot(Acond/max(Acond),'*b');    
                title('test_geomreconstruction');
                
                legend2 = sprintf('progressive rmsds (max=%g)',max(rmsds)); 
                %legend1 = sprintf('rec length deviations (max=%g)',max(abs(lengthdeviations)));                
                legend3 = sprintf('rec matrix conditions (max=%g)',max(Acond));
                
                legend({legend2,legend3});
                hold off;
            end
        end
        
        function [passed, rmsd1, rmsd2] = test_roundtrip(detailed)            
            dprintf('- Coords.test_roundtrip');
            
            if (nargin<1), detailed = 0; end;
            
            for i=1:data.GeomDB.nch,
                ch = data.GeomDB.ch(i);            
                if (mod(i,100) == 0 && detailed),
                    disp(sprintf('testing chain %i',i));
                end
                
                % bond -> geometry -> bond -> coords
                reccoords = ch.geometry.getcoords;                
                rmsd1(i) = reccoords.registerto(ch);
                
                % result1 -> ztransform -> betatransform -> geometry ->
                % bond -> coords
                reccoords.geometry.ztransform(protein.HEFG.HEFGSClassifier(1),1);
                reccoords.geometry.betatransform(protein.HEFG.HEFGSClassifier(1),1);
                
                cp = reccoords.geometry.getdata('cp');
                c = [reccoords.geometry.getdata('c');cp(end)];
                beta = reccoords.geometry.getdata({'betafromz'});
                
                transbond = geom.geometry2bond(c,cos(beta),sin(beta),reccoords.geometry.len);
                transcoords = geom.Coords(transbond);
                
                rmsd2(i) = transcoords.registerto(ch);
            end
            
            if (detailed)
                hist([rmsd1 rmsd2]);
            end
            
            passed = max(rmsd1 < 1e-5) && max(rmsd2 < 1e-5);
            
            if passed,    
                dprintf('PASSED.');
            else
                dprintf('PROBLEMS FOUND!');
            end
        end
    end
    
end
