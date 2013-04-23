%flips image along its vertical axis
%enter mirrorimage = mirror(pre)
function mimage = mirrory(coords)

sizec = size(coords,1);
%find mirror image of pre
mimage = zeros(sizec,3);
%find furthest aa on left and on right, to find midline
yleft = coords(1,2);
yright = coords(sizec,2);
for i = 1:sizec
    if coords(i,2) < yleft
        yleft = coords(i,2);
    end
    if coords(i,2) > yright
        yright = coords(i,2);
    end
end
mline = (yright + yleft)/2;
%flip all points to other side of midline
for i = 1:sizec
    y = 2*mline - coords(i,2);
    mimage(i,:) = [coords(i,1),y,coords(i,3)];
end
