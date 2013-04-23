function [] = ribbon3dnice( coords, cl, data, embedh )
%RIBBON3DNICE( coords, cl, data )
%   coords(i,:)         ith c-alpha
%   cl(:,i) (optional)  i-th fragment class set
%   data(optional)      struct with protein info
%                           .seq    sequence
%                           .name   name
%                           .res    resolution
%                           .bond   bonds
n = size(coords,1);

if (nargin < 2 || size(cl,1) == 0)
    cl = ones(size(coords,1)-3,1);
end

if (nargin < 3),
    data = [];
end

bond = geom.coords2bond(coords);
avgbond = mean(sqrt(sum(bond.^2,2)));

%% options
figure;
clf;
hold on;
ribbonaxh = gca;

width = 0.3;        % ribbon width, in units of average bond length
npoints = 10;       % number of interpolation points
facecolor = 'red';  % color of ribbon (if not class coloured)
facealpha = 0.85;   % alpha of ribbon
edgecolor = 'black';% color of ribbon edges
edgealpha = 0.15;   % alpha -.-

sphgrid = npoints;       % spheres resolution (number of subdivisions)
spheresalpha = 0.3;  % spheres alpha
spherescolor = 'red';
spheressize = width * avgbond;

width = width * avgbond;

%% create neighbor splines
centers = (coords(1:end-2,:) + coords(3:end,:))/2; % (n-2)x3

normv = zeros(n-2,3); % normal vector
zero = [];
for i = 1:n-2,
    vn = coords(i+1,:)-centers(i,:);
    if (norm(vn)>1e-10),
        normv(i,:) = cross(coords(i+2,:)-coords(i,:),vn); % (c-a)x(b-d)
        normv(i,:) = normv(i,:)/norm(normv(i,:));   
    else
        zero = [zero i];
    end
end
% check if we have any zeros
if (length(zero) == n-2)
    % choose arbitrary normal direction
    vn = rand(1,3);
    v = cross(coords(2,:)-coords(1,:),vn);
    normv = repmat(v/norm(v),n-2,1);
else
    for i=zero,
        %TODO: interpolate from next known directions
    end
end
% set first and last same as second and second last
normv = [normv(1,:); normv; normv(end,:)];

pos = coords + width*normv;
neg = coords - width*normv;

% create splines
pppos = spline(1:n, pos');
interpppos = ppval(pppos,1:1/npoints:n)';

ppneg = spline(1:n, neg');
interppneg = ppval(ppneg,1:1/npoints:n)';

pp = spline(1:n, coords');
interpp = ppval(pp,1:1/npoints:n)';

%% plot ribbon
%vertices
vertices = zeros(size(interppneg,1)+size(interpppos,1),3);
j = 1;
for i=1:2:size(vertices,1),
    vertices(i,:) = interpppos(j,:);
    vertices(i+1,:) = interppneg(j,:);
    j = j + 1;
end;

% faces
faces = zeros(size(vertices,1)/2-1,4);
for i=1:size(faces,1),
    faces(i,:) = [1,2,4,3] + (i-1)*2;
end

ribbonh = patch('Vertices',vertices,'Faces',faces,'FaceAlpha',facealpha,'EdgeColor',edgecolor,'EdgeAlpha',edgealpha);

userdata.currentcls = 1;
set(ribbonh,'UserData',userdata);

%% set colors

[colors, map] = setcolors(1);

%% plot center of mass
[X,Y,Z] = sphere(sphgrid);
com = sum(coords)/n;
surf(X*0.5+com(1),Y*0.5+com(2),Z*0.5+com(3),'FaceColor','black','EdgeColor','none', 'FaceAlpha',0.8);

setcameratarget(com,gca);

%% some more graphics options
set(gcf,'DoubleBuffer','on');
%set(gcf,'Color',[0,0,0]);
%set(gca,'Projection','perspective');

axis off;
axis equal;

lighth = camlight;
lighting gouraud;
material shiny;

set(gca,'CameraViewAngleMode','manual');
va = get(gca,'CameraViewAngle');
set(gca,'CameraViewAngle',0.8*va);
camh = cameratoolbar;
cameratoolbar('SetCoordSys','none') 

%set(gca,'Clipping','off');
%set(ribbonh,'AmbientStrength',0.6);
%set(ribbonh,'SpecularStrength',2);
set(findobj(gca,'-property','DiffuseStrength'),'DiffuseStrength',1);
set(findobj(gca,'-property','SpecularExponent'),'SpecularExponent',15);
set(findobj(gca,'-property','SpecularColorReflectance'),'SpecularColorReflectance',0.8);
set(findobj(gca,'-property','BackFaceLighting'),'BackFaceLighting','lit');
%set(gcf,'ButtonDownFcn',@setRotationCenter);
%set(gcf, 'WindowButtonDownFcn', @setRotationCenter);
%set(findobj('-property','WindowButtonDownFcn '),'WindowButtonDownFcn ',@setRotationCenter);

% callbacks
set(gcf, 'windowbuttondownfcn', {@btndownfcn, gca, [coords;com]})
set(gcf, 'windowbuttonupfcn', {@btnupfcn}); 
set(gcf, 'windowbuttonmotionfcn', {@btnmotionfcn}); 
set(gcf,'KeyPressFcn',{@keypressfcn, gca, [coords;com]});
vis.SetMouseWheel(@mouswheelfcn,gcf);

drawnow;
hold off;

%setup frag view axes
% fragaxh = axes('Position',get(gca,'Position'));
% set(fragaxh,'OuterPosition',[-0.1,0.6,0.5,0.5]);
% set(fragaxh,'HitTest','off');
% set(fragaxh,'Clipping','off');
% %set(fragaxh,'CameraTarget',[0,0,0]);
% %set(fragaxh,'CameraUpVector',[0,0,1]);
% %set(fragaxh,'CameraPosition',[-20,-20,5]);
% %set(fragaxh,'CameraTargetMode','auto');
% %set(fragaxh,'CameraPositionMode','auto');
% %set(fragaxh,'CameraUpvectorMode','auto');
% %set(fragaxh,'CameraViewAngleMode','auto');
% %set(fragaxh,'DataAspectRatio',[-1,1,1]); % swap x and y axis
% set(fragaxh,'View',[120,15]);
% camlight('right');
% lighting gouraud;
% material shiny;
% axis off;
% set(fragaxh,'Projection','Perspective');
% set(fragaxh,'DataAspectRatio',[1,1,1]);
% set(fragaxh,'XDir','reverse');
% set(fragaxh,'YDir','reverse');

% fragment view planes
% planev = [[0,0,0];[1,0,0];[1,1,0];[0,1,0];[0,0,1];[1,0,1];[0,1,1]] * avgbond;
% planef = [[1,2,3,4];[1,2,6,5];[1,4,7,5]];
% fragaxisplaneh = patch('Vertices',planev,'Faces',planef,'Tag',...
%     'fragaxplane','FaceAlpha',0.2,'Clipping','off','HitTest','off','FaceColor','White','Visible','off');


% fragment info axes
% fraginfoph = uipanel;                       
% set(fraginfoph,'Position',[0.8,0.9,0.2,0.1],'BackgroundColor',[0,0,0]);
% fraginfoaxh = axes('parent',fraginfoph);
% set(fraginfoaxh,'YDir','reverse');
% set(fraginfoaxh,'Visible','off');
% set(fraginfoaxh,'Tag','fraginfoaxh');
% axes(fraginfoaxh);
% fraginfotexth = text(-0.1,0.45,'no amino acid selected.','Color',[1,1,1],'Tag','fraginfotext','HorizontalAlignment','left');

% protein info axes
% protinfoph = uipanel;                       
% set(protinfoph,'Position',[0.8,0,0.2,0.1],'BackgroundColor',[0,0,0]);
% protinfoaxh = axes('parent',protinfoph);
% set(protinfoaxh,'YDir','reverse');
% set(protinfoaxh,'Visible','off');
% set(protinfoaxh,'Tag','protinfoaxh');
% axes(protinfoaxh);
% protinfotexth = text(-0.1,0.45,getprotinfo,'Color',[1,1,1],'Tag','protinfotext','HorizontalAlignment','left');            

%% toolbar
% toolbarh = uipanel;                       
% set(toolbarh,'Position',[0.95,0.12,0.05,0.1],'BackgroundColor',[0,0,0]);
% colorbtnh = uicontrol(... 
%                    'Style','togglebutton',...
%                    'Parent', toolbarh, ...
%                    'Units','normalized',...
%                    'HandleVisibility','callback', ...
%                    'Position',[0 0 1 0.2],...
%                    'String','red',...
%                    'Callback', @colorbtnfcn);
% cordtoggleh = uicontrol(... 
%                    'Style','togglebutton',...
%                    'Parent', toolbarh, ...
%                    'Units','normalized',...
%                    'HandleVisibility','callback', ...
%                    'Position',[0 0.25 1 0.2],...
%                    'String','cord',...
%                    'Callback', @cordtoggle);
% cordtoggleh = uicontrol(... 
%                    'Style','togglebutton',...
%                    'Parent', toolbarh, ...
%                    'Units','normalized',...
%                    'HandleVisibility','callback', ...
%                    'Position',[0 0.5 1 0.2],...
%                    'String','labels',...
%                    'Callback', @labelstoggle);
%                
% cordtoggleh = uicontrol(... 
%                    'Style','togglebutton',...
%                    'Parent', toolbarh, ...
%                    'Units','normalized',...
%                    'HandleVisibility','callback', ...
%                    'Position',[0 0.75 1 0.2],...
%                    'String','calphas',...
%                    'Callback', @spherestoggle);
               
%% class set switches
% if (size(cl,2) > 1)               
%     nextcls = uicontrol(...    
%                        'Parent', gcf, ...
%                        'Units','normalized',...
%                        'HandleVisibility','callback', ...
%                        'Position',[0.765 0.001 0.025 0.033],...
%                        'String','>',...
%                        'Callback', @nextclsfcn,'Enable','on');
%     prevcls = uicontrol(...
%                        'Parent', gcf, ...
%                        'Units','normalized',...
%                        'HandleVisibility','callback', ...
%                        'Position',[0.74 0.001 0.025 0.033],...
%                        'String','<',...
%                        'Callback', @prevclsfcn,'Enable','off','HitTest','off');
% end               
%set(cbarh,'OuterPosition',[0.1,0.001,0.7,0.033]);

axes(ribbonaxh); % set axes back to ribbon


%% callbacks
    function btndownfcn(src,eventdata, fh, coords)

        type = get(gcf,'SelectionType'); % normal = left, extend = middle, alt = right

        if (strcmp(type,'extend'))
            p = get(fh,'currentpoint');
            clickpoint = (p(1,:)+p(2,:))/2;        

            setcenteraa(getnext(coords,clickpoint),coords,fh);
        else
            cameratoolbar('down');            
        end
    end

    function btnupfcn(src,eventdata)
       type = get(gcf,'SelectionType');
       set(protinfotexth,'String',getprotinfo);
        
       if (strcmp(type,'extend'))

       else
            cameratoolbar('up')
       end
    end

    function btnmotionfcn(src, evt)
        if (exist('protinfotexth')) % startup
            set(protinfotexth,'String',getprotinfo);
        end
    end

    function str1 = getprotinfo
        [az,el] = view(ribbonaxh);
        i = 1;
        if (isstruct(data))
            str1(1) = {sprintf('Protein: %s',data.name)};            
            str1(2) = {sprintf('Length: %i',length(data.seq))};
            str1(3) = {sprintf('Resolution: %f',data.res)};
            str1(4) = {sprintf('--------------------------')};
            i = 5;
        end

        str1(i) = {sprintf('Az/El:  %2.0f / %2.0f',az,el)}; 
    end

    function keypressfcn(src, evnt, axh, coords)
    %     if ~isempty(evnt.Modifier)
    %          for ii = 1:length(evnt.Modifier)
    %             out = sprintf('Character: %c\nModifier: %s\nKey: %s\n',evnt.Character,evnt.Modifier{ii},evnt.Key);
    %             disp(out)
    %          end
    %       else
    %          out = sprintf('Character: %c\nModifier: %s\nKey: %s\n',evnt.Character,'No modifier key',evnt.Key);
    %          disp(out)
    %       end
        n = size(coords,1);
        current = camtarget(axh);
        if strcmp(evnt.Key,'leftarrow'),
            i = getnext(coords,current);        
            if (i == 1),
                i = n; 
            else
                i = i -1;
            end
            setcenteraa(i,coords,axh);
        else if strcmp(evnt.Key,'rightarrow'),
                i = getnext(coords,current);        
                if i == n, 
                    i = 1; 
                else
                    i = i + 1;
                end
                setcenteraa(i,coords,axh);
            else if strcmp(evnt.Key,'home'),
                setcenteraa(1,coords,axh);
                else
                   cameratoolbar('keypress');
                end
            end
        end

    end

    function mouswheelfcn(src,eventdata)
        scr = get(eventdata,'UnitsToScroll');
        camzoom(1-scr/9)
    end

    function colorbtnfcn(src,evt)

        val = get(src,'Value');
        if (val)
            set(ribbonh,'FaceColor',facecolor);
        else
            set(ribbonh,'FaceColor','interp');
        end

    end

    function prevclsfcn(src,evt)
        u = get(ribbonh,'userdata');
        c = u.currentcls;
        
        if c>1,
            colors = setcolors(c-1);
            u.currentcls = c -1;
            set(ribbonh,'UserData',u);
            set(nextcls,'Enable','on','HitTest','on');
        end
        if (c-1 == 1)
            set(prevcls,'Enable','off','HitTest','off');            
        end
    end

    function nextclsfcn(src,evt)
        u = get(ribbonh,'userdata');
        c = u.currentcls;
        
        if c<size(cl,2),
            colors = setcolors(c+1);
            u.currentcls = c+1;
            set(ribbonh,'UserData',u);
            set(prevcls,'Enable','on','HitTest','on');                        
        end        
        if (c+1 == size(cl,2))
            set(nextcls,'Enable','off','HitTest','off');            
        end
    end

%% frament axis
    function drawfrag(fragnum)        
        
        deleteall(fragaxh,'fragaxribbon');
        deleteall(fragaxh,'fragaxsphere');
        
        n = size(coords,1); %naa
        if (fragnum > 0 && fragnum < n-2)
            set(fragaxisplaneh,'Visible','on');
            fverts = fragverts(fragnum+1)';
            fverts = vertices(fverts,:);
            
            fragcolor = colors(2*fragnum*npoints+1,:);
            % normalize fragment rotation and transalte
            bond = coords(fragnum+1:fragnum+3,:)-coords(fragnum:fragnum+2,:);
            [bondr, rm] = geom.bond_normalizerotation(bond);
            faa = coords(fragnum,:); % first aa
            for i=1:size(fverts,1),
                fverts(i,:) = (fverts(i,:) - faa)*rm;
            end                        
            % rotate aa centers
            aacenters = coords(fragnum+[0:3],:);
            for i=1:size(aacenters,1),
                aacenters(i,:)=(aacenters(i,:) - faa)*rm;
            end
            
            % faces
            ffaces = zeros(size(fverts,1)/2-1,4);
            for i=1:size(ffaces,1),
                ffaces(i,:) = [1,2,4,3] + (i-1)*2;
            end
            
            hold on;
            axes(fragaxh);
            fragribbonh = patch('Vertices',fverts,'Faces',ffaces,'Tag','fragaxribbon','FaceColor',fragcolor,...
                'Clipping','off','HitTest','on','EdgeAlpha',edgealpha,'FaceAlpha',facealpha,'FaceLighting','gouraud');

            set(fragribbonh,'ButtonDownFcn',@fragBtnDwnFcn);
            
            % bonds
            for i=1:3,
                p = aacenters(i:i+1,:);
                line(p(:,1),p(:,2),p(:,3),'HitTest','off','Color',fragcolor,'LineWidth',3,'Tag','fragaxribbon');
            end
            %set(fragribbonh,'ButtonUpFcn',@fragBtnUpFcn);
            hold off;
            
            for j=1:4,
                drawsphere(aacenters(j,:),0.1,fragcolor,'fragaxsphere',fragaxh);
            end                       
                        
            % info text update
            str1(1) = {sprintf('Fragment No:  %i',fragnum)};
            str1(2) = {sprintf('Class:  %i',cl(fragnum))};
            if (isstruct(data)),
                str1(3) = {sprintf('Sequence: %i - %i - %i - %i',data.seq(fragnum),data.seq(fragnum+1),data.seq(fragnum+2),data.seq(fragnum+3))};
                [n1,l1,p1] = aaname3(data.seq(fragnum));
                [n2,l2,p2] = aaname3(data.seq(fragnum+1));
                [n3,l3,p3] = aaname3(data.seq(fragnum+2));
                [n4,l4,p4] = aaname3(data.seq(fragnum+3));
                str1(4) = {sprintf('Sequence: %s - %s - %s - %s',l1,l2,l3,l4)};
                str1(5) = {sprintf('Polarity: %s - %s - %s - %s',p1,p2,p3,p4)};
            end
            set(fraginfotexth,'String',str1);
            
            % set back axes
            axes(ribbonaxh);
        else % no fragment -> hide axis            
            set(fragaxisplaneh,'Visible','off');
        end
    end

    function fragBtnDwnFcn(scr, evt)
        disp('on');
    end

    function fragBtnUpFcn(scr, evt)
        disp('on');
    end

%% ribbon axis
    function setcenteraa(i, centers, axh) % centers and highlights frag(i-1,..,i+2)

        n = size(centers,1); % last is com        
        highlightcolors = colors; % get original colors

        % delete fragment aa spehres
        deleteall(axh,'fragsphere');
        
        centercolor = 'red';

        if (i < n-2 && i>1) % aa, light up fragment        

            fragcolor = colors(2*(i-1)*npoints+1,:);
            centercolor = fragcolor;
            
            fragidx = fragverts(i);
            highlightcolors(fragidx,:) = repmat(fragcolor,length(fragidx),1);                      

            highlightvert = fragidx([1:3:end]);
            highlightcolors(highlightvert,:) = repmat([1 1 1],length(highlightvert),1); % highligh frag

            % draw fragment aa spehres
            for j=[-1,1,2],
                drawsphere(centers(i+j,:),0.1,fragcolor,'fragsphere',axh);
            end
            
            % frag bonds
            axes(axh);
            for j=-1:1,
                p = centers(i+j:i+j+1,:);
                line(p(:,1),p(:,2),p(:,3),'HitTest','off','Color',fragcolor,'LineWidth',3,'Tag','fragsphere');
            end
            
        else if (i == n)
                set(fraginfotexth,'String','center of mass.');
            else if (i == n-1)
                set(fraginfotexth,'String','last amino acid.');
                else if (i == n-2)
                        set(fraginfotexth,'String','second last amino acid.');
                    else
                        set(fraginfotexth,'String','first amino acid.');
                    end
                end
            end
        end

        set(ribbonh,'FaceVertexCData',highlightcolors);
        setcameratarget(centers(i,:),axh, centercolor);
        drawfrag(i-1);
        axes(ribbonaxh);     
    end

    function setcameratarget(center, axh, color) % sets camera target and draws focus sphere

            if (nargin < 3),
                color = [];
            end

            % set camera target
            set(axh,'CameraTarget',center);

            %delete old sphere
            deleteall(axh,'centersphere');

            % draw sphere
            drawsphere(center, width, color, 'centersphere',axh);       
            axes(ribbonaxh);     
    end


%% helpers
    function [colors, map] = setcolors(index)
        
        ncl = max(max(cl)) + 1;
        map = [vis.class_colormap(ncl-1);[0,0,0]];
        colormap(map);
        
        colors = zeros(size(vertices,1),3) + max(cl(:,index))+1; %% no class = last in map = black
        for i=1:size(colors,1),
            if (i <= 2*npoints || i > (n-2)*2*npoints-2)
                continue
            end
            fragnum = floor((i+1)/(2*npoints));
            colors(i,:) = map(cl(fragnum,index),:);
        end
        
        set(ribbonh,'FaceColor','interp','FaceVertexCData',colors);
        
        % label colors
        if ~isempty(findall(ribbonaxh,'Tag','label')),
            deleteall(ribbonaxh,'label');
            labelstoggle([],[],index);
        end

%         % sphere colors
%         if (~isempty(findall(ribbonaxh,'tag','spheres')))
%             deleteall(ribbonaxh,'spheres');
%             spherestoggle([],[],index);
%         end

        % colorbar
        deleteall(gcf,'colorbartag');
        deleteall(gcf,'colorbarlabel');
        
        cbarh = colorbar('location','southoutside');
        
        set(cbarh,'OuterPosition',[0.1,0.001,0.7,0.033]);
        set(cbarh,'XTick',[]);
        set(cbarh,'HitTest','off');
        set(cbarh,'Layer','bottom');
        set(cbarh,'Tag','colorbartag');
        set(cbarh,'Clipping','off');        
        %get(cbarh)
        a = get(cbarh,'Children');
        set(a,'AlphaData',0.9);
        %get(a)
        axes(cbarh);
        xlim = get(cbarh,'XLim');
        if size(cl,2) > 1,
            text((-0.6)/size(map,1)*xlim(2),0.3,int2str(index),'Color','White','FontSize',14,'FontWeight','Bold');
        end

        f = freqs(cl(:,index), ncl);
        
        for i=1:size(map,1),           
            th = text((i-0.5)/size(map,1)*xlim(2),0.2,int2str(i));
            fth = text((i-0.5)/size(map,1)*xlim(2),2.2,int2str(f(i)),'Color','yellow');  
            if (i == size(map,1)) % last color is always black -> white font
                set(th,'Color','White');
            end
            set(th,'Tag','colorbarlabel');
            set(fth,'Tag','colorbarlabel');
            %set(th,'Parent',cbarh);
        end
            
        axes(ribbonaxh);                
    end

    function deleteall(baseh,tag) % delete all objects in baseh tree
        h = findall(baseh,'Tag',tag);
        if (h ~= 0),delete(h); end 
    end

    function [i] = getnext(points, point) % get points index next to point
        n = size(points,1);
        dist = points - repmat(point,n,1);
        dist = sum(dist.^2,2);
        [m,i] = min(dist);
    end

    function verts = fragverts(i) % get vertices belonging to frag i, where i = 2nd aa in frag
        ind = 1:2*3*npoints+4;
        start = 2*(i-2)*npoints-2;
        if (start < 0), start = 0;, end; % indexing glitch at aa 2
        verts = start + ind;
    end

    function h=drawsphere(center,radius,color,tag,axh)
               
        [X,Y,Z] = sphere(20);
        hold on;
        h = surf(axh,X*radius+center(1),Y*radius+center(2),Z*radius+center(3),'EdgeColor','none','FaceLighting','gouraud','BackFaceLighting','lit');
        if (nargin >= 3 && ~isempty(color))
            set(h,'FaceColor',color);
        end
        if (nargin >= 4)
            set(h,'Tag',tag);
        end
        hold off;      
    end

%% option handlers
    function cordtoggle(src,evt)
       if (isempty(findall(gcf,'Tag','cord')))
           hold on;   
           axes(ribbonaxh);
           plot3(interpp(:,1),interpp(:,2), interpp(:,3),'LineWidth',1,'Color','black','Tag','cord'); 
           hold off;
       else
           deleteall(gcf,'cord');
       end
       axes(ribbonaxh);
    end

    function labelstoggle(src,evt,cls)
        if nargin < 3,
            currentcls = get(ribbonh,'UserData');
            currentcls = currentcls.currentcls;
        else
            currentcls = cls;
        end
        
        if (isempty(findall(ribbonaxh,'Tag','label')))      
            for i=1:n-3,
                p = (coords(i+1,:)  + coords(i+2,:))/2;
                h = text(p(1), p(2), p(3), num2str(cl(i,currentcls)),'Color',map(cl(i,currentcls),:));
                set(h,'Tag','label');
            end
        else
            deleteall(ribbonaxh,'label');
        end
        axes(ribbonaxh);
    end

    function spherestoggle(src,evt,cls)
        %         if nargin < 3,
        %             currentcls = get(ribbonh,'UserData');
        %             currentcls = currentcls.currentcls;
        %         else
        %             currentcls = cls;
        %         end
        if (isempty(findall(ribbonaxh,'tag','spheres')))
            [X,Y,Z] = sphere(sphgrid);
            hold on;
            axes(ribbonaxh);
            for i=1:n,
                surf(X*spheressize+coords(i,1),Y*spheressize+coords(i,2),Z*spheressize+coords(i,3),...
                    'FaceColor',spherescolor,...
                    'EdgeColor','none',...
                    'FaceAlpha',spheresalpha,...
                    'FaceLighting','gouraud',...
                    'BackFaceLighting','lit',...
                    'Tag','spheres');
                %             if (i>2 && i<n-2)
                %             	set(spheresh(i),'FaceColor',map(ceil((cl(i-1,index) + cl(i-2,index))/2),:));
                %             else
                %             	set(spheresh(i),'FaceColor',map(max(cl(:,index)),:));
                %             end
            end
            material shiny;
            hold off;
        else
            deleteall(ribbonaxh,'spheres');
        end
        axes(ribbonaxh);
    end


end