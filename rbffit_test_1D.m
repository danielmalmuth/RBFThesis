%RBFFIT Fits the thin plate spline radial basis function (RBF) to the data (x(i),y(i),f(i)).
%    RBFFIT(x,y,f) finds the coefficients of the thin-plate spline RBF that
%    interpolates the data (x(i),y(i),f(i)), where x and y are the independet
%    variables.  The thin-plate spline is defined as
%                          phi(r) = (r^2)*log(r)
%
%   Inputs:
%          x : x coordinates of the given nodes.
%          f : Scalar function values to interpolate at each x(i)
%
%   Outputs:
%        lam : the coefficients (weights) for the 2-D thin plate spline
%              interpolant.
%
%  See also: rbfval.
function lam = rbffit_test_1D(x,f,solver,schur,trunc,deg)

if ~exist('solver','var')
    % fifth parameter does not exist, so default it to something
    solver = 'Backslash';
end

if ~exist('schur','var')
    % fifth parameter does not exist, so default it to something
    schur = true;
end

if ~exist('trunc','var')
    % fifth parameter does not exist, so default it to something
    trunc = 0;
end

sz = size(x);

% Flatten x, y, and f;
x = x(:); f = f(:);
n = length(x);

if n ~= length(x) || n ~=length(f)
    error('The sizes of the input vectors must be the same');
end

%
% Compute the pairwise distances between all the nodes.
%

% Prepare the distance squared matrix.
A = zeros(n);

[xd1,xd2] = meshgrid(x);  % x coordinates
A = A + (xd1-xd2).^2;

clear xd1 xd2;

% Need to apply the TPS radial function to A.  However, the TPS is r^2*log(r)
% and will have numerical trouble when r=0.  It should be 0 when r=0, so we
% explicity make this the case.
id = 1:(n+1):n^2;    % indicies for the diagonal of A (i.e. where r=0).
A(id) = 1;           % make the distance 1 since log(1) = 0.
% Compute the TPS interpolation matrix.  Note that we are dealing with the
% distances squared.  The TPS kernel with distances squared reduces to
% 0.5*A*log(A)
A = 0.5*A.*log(A);

% Add on the extra polynomial conditions.
ev = ones(n,1);

if deg == 0
    B = ev;
else
    B = [ev x];
    for i = 2:deg
        B = [B x.^i];
    end
end

if schur
    lambda = ones(size(B,2),1);
    S = B'*(A\B);
    switch solver
        case 'Backslash'
            lambda = S\(B'*(A\f));
        case 'Pseudoinverse'
            lambda = pinv(S,trunc)*B'*(A\f);
        case 'SVD'
            [U,S,V] = svd(B,0);
            L = S > trunc*S(1,1);
            bsize = size(B,2);
            l = diag(L);
            if (isequal(ones(bsize,1),l) == 0)
                newdim = sum(l(:)==1);
                U = U(1:n,1:newdim);
                S = S(1:newdim,1:newdim);
                V = V(1:bsize,1:newdim);
                lambda = lsqminnorm(U'*(A\(U*S*V')),U'*(A\f));
            else
                lambda = (U'*(A\(U*S*V')))\(U'*(A\f));
            end
        case 'QR'
            [Q,R] = qr(B);
            L = abs(R) > trunc*abs(R(1,1));
            bsize = size(B,2);
            l = diag(L);
            if (isequal(ones(bsize,1),diag(L)) == 0)
                newdim = sum(l(:)==1);
                Q = Q(1:n,1:newdim);
                R = R(1:newdim,1:bsize);
                lambda = lsqminnorm((R'*Q'*(A\Q)*R),(R'*Q'*(GaussPivot(A,f))));
            else
                lambda = (R'*Q'*(A\Q)*R)\(R'*Q'*(GaussPivot(A,f)));
            end
            % Why not use backslash for full rank?
        case 'GECP'
            lambda = GaussPivot(S,(B'*(A\f)));
    end
    c = A\(f-B*lambda);
    lam = [c;lambda];
else
    bsize = size(B,2);
    A = [[A B];[B' zeros(bsize)]];
    asize = size(A,2);
    
    switch solver
        case 'Backslash'
            lam = A\[f;zeros(bsize,1)];
        case 'Pseudoinverse'
            lam = pinv(A)*[f;zeros(bsize,1)];
        case 'SVD'
            [U,S,V] = svd(A);
            L = S > trunc*S(1,1);
            l = diag(L);
            if (isequal(ones(asize,1),l) == 0)
                newdim = sum(l(:)==1);
                U = U(1:asize,1:newdim);
                S = S(1:newdim,1:newdim);
                V = V(1:asize,1:newdim);
                %     lambda = (U'*(A\U)*S*V')\(U'*(GaussPivot(A,f)));
                %     lambda = (U'*(A\U)*S*V')\(U'*(A\f));
                lam = lsqminnorm((U*S*V'),[f;zeros(bsize,1)]);
            else
                lam = (U*S*V')\([f;zeros(bsize,1)]);
            end
        case 'QR'
            [Q,R] = qr(A);
            L = R > trunc*R(1,1);
            l = diag(L);
            if (isequal(ones(asize,1),diag(L)) == 0)
                newdim = sum(l(:)==1);
                Q = Q(1:asize,1:newdim);
                R = R(1:newdim,1:asize);
                lam = lsqminnorm((Q*R),[f;zeros(bsize,1)]);
            else
                lam = (Q*R)\([f;zeros(bsize,1)]);
            end
            % Why not use backslash for full rank?
        case 'GECP'
            lam = GaussPivot(A,[f;zeros(bsize,1)]);
    end
end