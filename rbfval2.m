%RBFVAL2 Evaluates the 2D thin-plate spline (TPS) radial basis function interpolant 
%       specified by the coefficients lambda at the points ux,uy.
%
% s = RBFVAL2(lam,x,y,ux,uy) evaluates the 2D thin plate spline interpolant
% specified by the coefficents lam(j) (obtained from rbffit) and the node 
% locations (x(j),y(j)) at the points (ux(i),uy(i)).
%
%   Inputs: 
%
%        lam : The thin-plate expansion coefficients computed from rbffit.
%          x : x coordinates of the given nodes.
%          y : y coordinates of the given nodes.
%         ux : x coordinates of the evaluation locations.
%         uy : y coordinates of the evaluation locations.
%
%  See also: rbffit
function p = rbfval2(lam,x,y,ux,uy,deg)

sz = size(ux);

% Flatten x, y, lam, ux, uy;
x = x(:); y = y(:); lam = lam(:)'; ux = ux(:); uy = uy(:);
n = length(x);
m = length(ux);

% Find number of terms in polynomial basis when given degree
polysize = 0;
for i = 1:(deg+1)
    polysize = polysize + i;
end


% The following is not the most computationally efficient way to do the
% evaluation.  Instead it is designed to be efficient in terms of memory
% use.

% Prepare the distance squared matrix.
B = zeros(n+polysize,1);
p = zeros(m,1);
for j=1:m
   B(1:n,1) = (x-ux(j)).^2 + (y - uy(j)).^2;
   
   % Need to apply the TPS radial function to B.  However, the TPS is r^2*log(r)
   % and will have numerical trouble when r=0.  It should be 0 when r=0, so we
   % explicity make this the case.
   B = 0.5*B.*log(B+eps);
   
   % Polynomial part.
   
   if deg == 0
       polybasis = [1];
   else
       polybasis = [1;ux(j);uy(j)];
       for i = 2:deg
           polybasis = [polybasis;ux(j).*polybasis((end-i+1):end);uy(j).^i];
       end
   end
   
%    B(n+1:n+6,1) = [1;ux(j);uy(j);ux(j).^2;uy(j).^2;ux(j).*uy(j)];
   B(n+1:n+polysize,1) = polybasis;
   
   % Evaluate the RBF
   p(j,1) = lam*B;
end

p = reshape(p,sz);



