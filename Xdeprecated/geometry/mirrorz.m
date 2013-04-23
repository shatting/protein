%flips image along its z-axis
%enter mirrorimage = mirror(pre)
function mimage = mirrorz(coords)

sizec = size(coords,1);
%find mirror image of pre
mimage = zeros(sizec,3);
%find furthest aa on top and on bottom, to find midline
xbottom = coords(1,3);
xtop = coords(sizec,3);
for i = 1:sizec
    if coords(i,3) < xbottom
        xbottom = coords(i,3);
    end
    if coords(i,3) > xtop
        xtop = coords(i,3);
    end
end
mline = (xtop + xbottom)/2;
%flip all points to other side of midline
for i = 1:sizec
    z = 2*mline - coords(i,3);
    mimage(i,:) = [coords(i,1),coords(i,2),z];
end
