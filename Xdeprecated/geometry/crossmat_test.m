function crossmat_test()


for i=1:100
    a = rand(1,3);
    b = rand(1,3);    
    mycorssc = a*geom.crossmat(b)
    matlabcross = cross(a,b)
    a*mycorssc'
    b*mycorssc'
end