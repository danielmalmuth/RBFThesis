function mat = localinterpmat_1D(x,deg)
    [xi,xj] = meshgrid(x,x);
    rbf_part = abs((xj-xi)).^3;
    
    if deg == 0
        poly = ones(length(x),1);
    else
    poly = [ones(length(x),1) x];
    for i = 2:deg
        poly = [poly x.^i];
    end
    mat = [rbf_part poly;poly' zeros(deg+1)];
    mat = decomposition(mat);
end