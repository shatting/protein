%flips image along its vertical axis
%enter mirrorimage = mirror(pre)
function mimage = mirrorx(coords)

sizec = size(coords,1);
%find mirror image of pre
mimage = zeros(sizec,3);
%find furthest aa on left and on right, to find midline
xleft = coords(1,1);
xright = coords(sizec,1);
for i = 1:sizec
    if coords(i,1) < xleft
        xleft = coords(i,1);
    end
    if coords(i,1) > xright
        xright = coords(i,1);
    end
end
mline = (xright + xleft)/2;
%flip all points to other side of midline
for i = 1:sizec
    x = 2*mline - coords(i,1);
    mimage(i,:) = [x,coords(i,2),coords(i,3)];
end
