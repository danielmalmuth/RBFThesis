function [xp,yp] = diagproj(x,y,alpha)

n = length(x);
x = x(:); y = y(:);
xp = zeros(n,1);
yp = zeros(n,1);

for i = 1:n
    co = dot([x(i);y(i)],[sqrt(1/2);sqrt(1/2)])/norm([x(i) y(i)]);
    si = sign(y(i)-x(i))*sqrt(1-co^2);
    xp(i) = norm([x(i) y(i)])*sqrt(1/2)*(co-si*alpha);
    yp(i) = norm([x(i) y(i)])*sqrt(1/2)*(co+si*alpha);
end

end