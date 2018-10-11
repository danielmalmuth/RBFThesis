function vec = localrbfinterp_1D(lam,x_star,x,deg)
    r = abs(x_star-x);
    rbf = r.^3;
    if deg == 0
        poly = 1;
    else
    poly = [1;x_star];
    for i = 2:deg
        poly = [poly;x_star.^i];
    end
    funcs = [rbf;poly];
    vec = dot(lam,funcs);
end