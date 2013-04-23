function [ classes_t, polys ] = classes2d_manual(data)
%INITIAL_CLASSES Summary of this function goes here
%   Detailed explanation goes here

n = size(data,1);
figure;
maximize_fig;

scatterplot([data(:,1),data(:,2)],ones(n,1),'','',1000000, 0.1, 0);
axis square;
%xlims = [min(p(:,1)) max(p(:,1))];
%ylims = [min(p(:,2)) max(p(:,2))];
%set(gcf,'WindowButtonMotionFcn',@mousemove);

toolbarh = uipanel;                       
set(toolbarh,'Position',[0.85,0,0.15,1],'BackgroundColor',[0,0,0]);

axis([0 1 0 1]);

hold on
% Initially, the list of points is empty.
% Loop, picking up the points.
disp('Left mouse button picks points.')
disp('Right mouse button picks last point of current polygon.')
disp('Middle mouse button finishes.')

p = 0;
but = 1;

while but ~= 2,
    xy = [];
    n = 0;
    but = 1;
    p = p + 1;
    
    while but == 1
        [xi,yi,but] = ginput(1);
        n=size(xy,1) + 1;
        if (but~=2)
            polys{p}.circle{n} = plot(xi,yi,'bo','Tag','circle');
        end
        xy(n,1) = xi;
        xy(n,2) = yi;
        if (size(xy,1) > 1),
            polys{p}.hline{n} = line(xy(end-1:end,1),xy(end-1:end,2));
        end
    end
    if(size(xy,1) < 3),
        continue;
    end
    
    polys{p}.polyxy = 2*(xy-0.5); %scatterplot is in [0,1]x[0,1], so shift
    
    %idx = convhull(x,y);
    %x = x(idx);
    %y = y(idx);
%     x = [xy(:,1);xy(1,1)];
%     y = [xy(:,2);xy(1,2)];
    polys{p}.hline{n+1} = line(xy([1,end],1),xy([1,end],2));        
    
    drawbuttons;
end
hold off;

classes_t = zeros(size(geometry_means,1),1);
shift = 0;
for p = 1:length(polys),
    idx = idx_inpoly(polys{p}.polyxy);
    
    if (sum(idx) == 0)
        shift = shift + 1; %% skip empty classes
        continue;
    else
        classes_t(idx) = p - shift;
    end
    polys{p-shift} = polys{p-shift}.polyxy;
end

polys = polys(1:end-shift);

ncl = length(polys);
classes_t(classes_t==0) = ncl + 1;
map = class_colormap(ncl + 1);

h = findall(gcf,'Tag','circle');
for i=1:length(h),
    delete(h(i));
end

scatterplot([geometry_means(:,3),geometry_means(:,4)],classes_t,'','',1000000, 0.1, 0);
colormap(map);
colorbar;

    function plot_all(),
        
        scatterplot([geometry_means(:,3),geometry_means(:,4)],classes_t,'','',1000000, 0.1, 0);        
    end
%     function mousemove(src, evt)
%         p = get(src,'CurrentPoint');
%         if (n>0)
%             line([xi,p(1)],[yi,p(2)]);
%         end
%         %disp('------------');
%         %get(evt)
%     end
    function idx = idx_inpoly(poly)
        xv = [poly(:,1);poly(1,1)];
        yv = [poly(:,2);poly(1,2)];
        idx = inpolygon(geometry_means(:,3),geometry_means(:,4),xv,yv);
    end

    function deleteall(baseh,tag) % delete all objects in baseh tree
        h = findall(baseh,'Tag',tag);
        if (h ~= 0),delete(h); end 
    end

    function drawbuttons()
        deleteall(gcf,'delete_button');
        deleteall(gcf,'up_button');
        deleteall(gcf,'down_button');
        deleteall(gcf,'label');
        np = length(polys);
        for i=1:np,
            if (p<=10),
                f = 0.1;
            else                
                f = 1/np;
            end

            pos =   [0.05,  1-i*f,          0.7,    0.9*f];
            posu =  [0.75,  1-i*f+0.5*f,    0.2,    0.4*f];
            posd =  [0.75,  1-i*f,          0.2,    0.4*f];               

            c = sum(idx_inpoly(polys{i}.polyxy));
            uicontrol(...                    
                       'Parent', toolbarh, ...
                       'Units','normalized',...
                       'HandleVisibility','callback', ...
                       'Position',pos,...
                       'String',sprintf('%i [%i]',i,c),...
                       'UserData',i,...
                       'Callback', {@deleteclass,i},...
                       'Tag','delete_button');
            if (i>1),
                uicontrol(...                    
                           'Parent', toolbarh, ...
                           'Units','normalized',...
                           'HandleVisibility','callback', ...
                           'Position',posu,...
                           'String','up',...
                           'UserData',i,...
                           'Callback', {@switchclass,i-1,i},...
                           'Tag','up_button');
            end
            if (i<np),
                uicontrol(...                    
                           'Parent', toolbarh, ...
                           'Units','normalized',...
                           'HandleVisibility','callback', ...
                           'Position',posd,...
                           'String','down',...
                           'UserData',i,...
                           'Callback', {@switchclass,i,i+1},...
                           'Tag','down_button');
            end                
                   
            text('Position',polys{i}.polyxy(1,:)/2+0.5+[0.01,0],'Tag','label','String',int2str(i),'BackgroundColor','white','FontSize',15);
        end
    end

    function switchclass(src,evt,a,b)
        if (length(xy) > 0)
            disp('please close current polygon before moving classes.');
            return;
        end
        t = polys{a};
        polys{a} = polys{b};
        polys{b} = t;
        drawbuttons;
    end

    function deleteclass(src,evt,np)
        
        if (length(xy) > 0)
            disp('please close current polygon before deleting.');
            return;
        end
        
        for n=1:length(polys{np}.hline)
            delete(polys{np}.hline{n});
        end
        for n=1:length(polys{np}.circle)
            delete(polys{np}.circle{n});
        end
        
        for i=np:length(polys)-1
            polys{i} = polys{i+1};
        end
        polys = polys(1:end-1);
        
        p = length(polys) + 1;
        
        drawbuttons;
    end

end